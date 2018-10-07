#ifndef SU2_REALVECSLICE_HPP
#define SU2_REALVECSLICE_HPP

#include "exprT.hh"

namespace Common {

  template  <typename T, int N> class RealVecSlice;

  template <typename T, int N> std::ostream& operator<<(std::ostream& out, const RealVecSlice<T,N>& v);
  template <typename T, int N> std::istream& operator>>(std::istream& in, RealVecSlice<T,N>& v);


/*
 * Definition of a class RealVecSlice that implements an expression template technique
 * with statically defined size
 * @author Andrea Lani
 */
template <typename T, int N = 0>
class RealVecSlice : public ETVEC(ExprT,RealVecSlice,T,N) {
public:

  /// Default constructor
   RealVecSlice(T* data) : ETVEC(ExprT,RealVecSlice,T,N)(this), m_data(data) {}

  /// Copy constructor from other slice
   RealVecSlice(const RealVecSlice<T,N>& in) : ETVEC(ExprT,RealVecSlice,T,N)(this), m_data(in.m_data) {}

  /// Default destructor
   ~RealVecSlice() {}

  /// Overloading of assignment operator(s)
#define VECSLICE_EASSIGN_OP(__op__) \
  template <typename EXPR>						\
   const RealVecSlice<T,N>& operator __op__ (EETYPEV(EXPR) expr) \
    {EVECLOOP(i, 0, N, m_data[i] __op__ expr.at(i)); return *this;}
VECSLICE_EASSIGN_OP(=)
VECSLICE_EASSIGN_OP(+=)
VECSLICE_EASSIGN_OP(-=)
VECSLICE_EASSIGN_OP(*=)
VECSLICE_EASSIGN_OP(/=)
#undef VECSLICE_EASSIGN_OP

  /// Overloading of assignment operator(s) with constants
#define VECSLICE_EASSIGN_OP_CONST(__op__) \
   const RealVecSlice<T,N>& operator __op__ (T expr) \
  {EVECLOOP(i, 0, N, m_data[i] __op__ expr); return *this;}
VECSLICE_EASSIGN_OP_CONST(=)
VECSLICE_EASSIGN_OP_CONST(+=)
VECSLICE_EASSIGN_OP_CONST(-=)
VECSLICE_EASSIGN_OP_CONST(*=)
VECSLICE_EASSIGN_OP_CONST(/=)
#undef VECSLICE_EASSIGN_OP_CONST

  /// Overloading of the assignment operator "=" with RealVecSlice
   const RealVecSlice<T,N>& operator= (const RealVecSlice<T,N>& other)
  {
    cf_assert(&other != this);
    for (std::size_t i = 0; i < N; ++i) {
      m_data[i] = other.m_data[i];
    }
    return *this;
  }

  /// return the array size
   std::size_t size() const {return N;}

  /// copy content of another array
  template <typename ARRAY>
   void copyFrom(const ARRAY& in){EVECLOOP(i, 0, N, m_data[i] = in[i]);}

  /// copy content of another array
  template <typename ARRAY>
   void copyTo(ARRAY& in) {EVECLOOP(i, 0, N, in[i] = m_data[i]);}

  /// @return the raw data
   T* ptr() {return m_data;}

  /// Overloading of the stream operator "<<" for the output.
  /// "\n"ine introduced at the end of every line of the matrix.
  friend std::ostream& operator<<
    LTGT (std::ostream& out, const RealVecSlice<T,N>& v);

  /// Overloading of the stream operator ">>" for the input
  friend std::istream& operator>>
    LTGT (std::istream& in, RealVecSlice<T,N>& v);

private:

  /// array data
  T* m_data;
};

//////////////////////////////////////////////////////////////////////////////

/// Definition of a class RealVecSlice that implements an expression template technique
/// with dynamical size
/// @author Andrea Lani
template <typename T>
class RealVecSlice<T,0> : public ETVEC(ExprT,RealVecSlice,T,0) {
public:

  // Constructor from preallocated memory
   RealVecSlice(T* data, std::size_t ns) :
    ETVEC(ExprT,RealVecSlice,T,0)(this), m_size(ns), m_data(data) {}

  /// Copy constructor from other slice
   RealVecSlice(const RealVecSlice<T,0>& in) :
    ETVEC(ExprT,RealVecSlice,T,0)(this), m_size(in.size()), m_data(in.m_data) {}

  /// Default destructor
   ~RealVecSlice() {}

   /// Overloading of assignment operator(s)
#define VECSLICE0_EASSIGN_OP(__op__) \
  template <typename EXPR>						\
   const RealVecSlice<T,0>& operator __op__ (EETYPEV(EXPR) expr) \
    {EVECLOOP(i, 0, EGETSIZE1(ENMAX(0,EXPR::SIZE1)), m_data[i] __op__ expr.at(i)); return *this;}
VECSLICE0_EASSIGN_OP(=)
VECSLICE0_EASSIGN_OP(+=)
VECSLICE0_EASSIGN_OP(-=)
VECSLICE0_EASSIGN_OP(*=)
VECSLICE0_EASSIGN_OP(/=)
#undef VECSLICE0_EASSIGN_OP

/// Overloading of assignment operator(s) with constants
#define VECSLICE0_EASSIGN_OP_CONST(__op__) \
   const RealVecSlice<T,0>& operator __op__ (T expr) \
  {EVECLOOP(i, 0, size(), m_data[i] __op__ expr); return *this;}
VECSLICE0_EASSIGN_OP_CONST(=)
VECSLICE0_EASSIGN_OP_CONST(+=)
VECSLICE0_EASSIGN_OP_CONST(-=)
VECSLICE0_EASSIGN_OP_CONST(*=)
VECSLICE0_EASSIGN_OP_CONST(/=)
#undef VECSLICE0_EASSIGN_OP_CONST

  /// Overloading of the assignment operator "=" with RealVecSlice
   const RealVecSlice<T,0>& operator= (const RealVecSlice<T,0>& other)
  {
    cf_assert(&other != this);
    for (std::size_t i = 0; i < m_size; ++i) {
      m_data[i] = other.m_data[i];
    }
    return *this;
  }

  /// return the array size
   std::size_t size() const {return m_size;}

  /// @return the raw data
   T* ptr() {return m_data;}

  /// Overloading of the stream operator "<<" for the output.
  /// "\n"ine introduced at the end of every line of the matrix.
  friend std::ostream& operator<<
    LTGT (std::ostream& out, const RealVecSlice<T,0>& v);

  /// Overloading of the stream operator ">>" for the input
  friend std::istream& operator>>
    LTGT (std::istream& in, RealVecSlice<T,0>& v);

private:

  /// array size
  std::size_t m_size;

  /// array data
  T* m_data;
};

//////////////////////////////////////////////////////////////////////////////

template <typename T, int N>
std::ostream& operator<< (std::ostream& out, const RealVecSlice<T,N>& v)
{
  const std::size_t size = v.size();
  for (std::size_t i = 0; i < size; ++i)
    out << v.m_data[i] << " " ;
  return out;
}

//////////////////////////////////////////////////////////////////////////////

template <typename T, int N>
std::istream& operator>> (std::istream& in, RealVecSlice<T,N>& v)
{
  const std::size_t size = v.size();
  for (std::size_t i = 0; i < size; ++i)
    in >> v.m_data[i];
  return in;
}


} /*--- End of namespace Common ---*/


#endif
