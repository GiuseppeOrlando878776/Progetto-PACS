#ifndef SU2_REALMATRIX_HPP
#define SU2_REALMATRIX_HPP

//#include "MathTools/CFMatSlice.hh"
#include "MathTools/RealVec.hh"


namespace Common {

  template  <typename T, int N, int M> class RealMatrix;
  template  <typename T, int N> class RealVec;

  template <typename T, int N, int M> std::ostream& operator<<(std::ostream& out, const RealMatrix<T,N,M>& A);
  template <typename T, int N, int M> std::istream& operator>>(std::istream& in, RealMatrix<T,N,M>& A);

/*
  * Definition of a class RealMatrix that implements an expression template technique
  * with statically defined size
  * @author G. Orlando
  */
template <typename T, int N = 0, int M = 0>
class RealMatrix : public ETMAT(MatExprT,RealMatrix,T,N,M) {
public:

  /// Default constructor
   RealMatrix() : ETMAT(MatExprT,RealMatrix,T,N,M)(this) {this->initialize(T());}

  /// Copy constructor from expression
  template <typename EXPR>
   RealMatrix(MEETYPEV(EXPR) expr) : ETMAT(MatExprT,RealMatrix,T,N,M)(this)
  {EMATLOOP(i,0,N,j,0,M, (*this)(i COMMA j) = expr.at(i COMMA j));}

  /// Copy constructor from RealMatrix
   RealMatrix(const RealMatrix<T,N,M>& orig) : ETMAT(MatExprT,RealMatrix,T,N,M)(this)
  {EMATLOOP(i,0,N,j,0,M, (*this)(i COMMA j) = orig(i COMMA j));}

  /// Default destructor
   ~RealMatrix() {}

  /// Overloading of assignment operator(s)
#define MAT_EASSIGN_OP(__op__) \
  template <typename EXPR>						\
     const RealMatrix<T,N,M>& operator __op__ (MEETYPEV(EXPR) expr) \
    {EMATLOOP(i,0,N,j,0,M, (*this)(i COMMA j) __op__ expr.at(i COMMA j)); return *this;}
MAT_EASSIGN_OP(=)
MAT_EASSIGN_OP(+=)
MAT_EASSIGN_OP(-=)
MAT_EASSIGN_OP(*=)
MAT_EASSIGN_OP(/=)
#undef MAT_EASSIGN_OP

   /// Overloading of assignment operator(s) with constants
#define MAT_EASSIGN_OP_CONST(__op__) \
   const RealMatrix<T,N,M>& operator __op__ (T expr) {EVECLOOP(i, 0, size(), m_data[i] __op__ expr); return *this;}
MAT_EASSIGN_OP_CONST(=)
MAT_EASSIGN_OP_CONST(+=)
MAT_EASSIGN_OP_CONST(-=)
MAT_EASSIGN_OP_CONST(*=)
MAT_EASSIGN_OP_CONST(/=)
#undef MAT_EASSIGN_OP_CONST

 /// Overloading of the assignment op=() operators for RealVec
#define MAT_EASSIGN_OP_VEC(__op__)  				\
  const RealMatrix<T,N,M>& operator __op__ (const RealVec<T,N>& diag)		\
{                                                            \
  for (std::size_t i = 0; i < N; ++i) { \
    (*this)(i,i) __op__ diag[i];     \
  }                                  \
  return *this;                      \
}
  MAT_EASSIGN_OP_VEC(=)
  MAT_EASSIGN_OP_VEC(+=)
  MAT_EASSIGN_OP_VEC(-=)
  MAT_EASSIGN_OP_VEC(*=)
  MAT_EASSIGN_OP_VEC(/=)
#undef MAT_EASSIGN_OP_VEC

  /// copy content of another array
  template <typename ARRAY>
   void copyFrom(const ARRAY& in) {EVECLOOP(i, 0, size(), m_data[i] = in[i]);}

  /// copy content of another array
  template <typename ARRAY>
   void copyTo(ARRAY& in) {EVECLOOP(i, 0, size(), in[i] = m_data[i]);}

  /// @return a vector slice with fixed size
  template <int NV, int MV>
   RealMatrixSlice<T,NV,MV> slice(std::size_t nn, std::size_t mm)
    {return RealMatrixSlice<T,NV,MV>(N, M, &m_data[nn*M+mm]);}

  /// @return a vector slice
   RealMatrixSlice<T,0,0> slice(std::size_t nn, std::size_t mm, std::size_t ns, std::size_t ms)
  {return RealMatrixSlice<T,0,0>(&m_data[nn*M+mm], N, M, ns, ms);}

  /// @return the raw data
   T* ptr() {return &m_data[0];}

  /// return the array size
   std::size_t size() const {return N*M;}

  /// return the number of rows
   std::size_t nbRows() const {return N;}

  /// return the number of columns
   std::size_t nbCols() const {return M;}

  /// Overloading of the operator"()" for assignment.
   T& operator() (std::size_t i, std::size_t j) {assert(i < N); assert(j < M); return m_data[i*M+j];}

  /// Overloading of the operator"()" (doesn't allow assignment)
   T operator() (std::size_t i, std::size_t j) const {assert(i < N); assert(j < M);return m_data[i*M+j];}

  /// Overloading of the stream operator "<<" for the output.
  /// "\n"ine introduced at the end of every line of the matrix.
  friend std::ostream& operator<<
    LTGT (std::ostream& out, const RealMatrix<T,N,M>& A);

  /// Overloading of the stream operator ">>" for the input
  friend std::istream& operator>>
    LTGT (std::istream& in, RealMatrix<T,N,M>& A);

  /// Accessor to individual entry
  /// @param i row index
  /// @param j column index
   T at (std::size_t i, std::size_t j) const {return (*this)(i,j);}

  /// Calculate the determinant of a 2*2 matrix
   T determ2() const;

  /// Calculate the determinant of a 3*3 matrix
   T determ3() const;

  /// Calculate the determinant of a 4*4 matrix
   T determ4() const;

  /// Invert a diagonal matrix
   void invertDiag(RealMatrix<T,N,M>& result) const;

  /// Returns the transposes of the object matrix.
  /// The object remains untransposed.
  /// @param result transposed matrix
   void transpose(RealMatrix<T,0>& result) const
  {EMATLOOP(j, 0, M, i, 0, N, result(i COMMA j) = m_data[j*M + i]);}

  /// Puts the copy of the row of the RealMatrix in a supplied vector
  /// @param iRow  number of the row
  template <typename VEC>
   void putRow(std::size_t iRow, VEC& v) const
  {
    assert(iRow < nbRows()); assert(v.size() == nbCols());
    EVECLOOP(i, 0, nbCols(), v[i] = m_data[iRow*nbCols() + i]);
  }

  /// Puts the copy of the column of the RealMatrix in a supplied vector
  /// @param iCol  number of the column
  template <typename VEC>
   void putColumn(std::size_t iCol, VEC& v) const
  {
    assert(iCol < nbCols()); assert(v.size() == nbRows());
    EVECLOOP(i, 0, nbRows(), v[i] = m_data[nbCols()*i + iCol]);
  }

  /// Set a row of the matrix.
  /// @param row   vector with the row to set
  /// @param iRow  number of the row
  template <typename VEC>
   void setRow(const VEC& row, std::size_t iRow)
  {
    assert(row.size() == nbCols());
    std::size_t istart = iRow*nbCols();
    for (std::size_t i = 0; i < nbCols(); ++i, ++istart) {
      m_data[istart] = row[i];
    }
  }

  /// Set a column of the matrix.
  /// @param col   vector with the column to set
  /// @param iCol  number of the column
  template <typename VEC>
   void setColumn(const VEC& col, std::size_t iCol)
  {
    assert(col.size() == nbRows());
    for (std::size_t i = 0; i < nbRows(); ++i) {
      m_data[iCol+nbCols()*i] = col[i];
    }
  }

  /// Gets a Vec with the copy of the row of the RealMatrix.
  /// @param iRow  number of the row
  /// @return the specified row of the matrix.
  template <typename VEC>
   VEC getRow(std::size_t iRow) const
  {assert(iRow < nbRows()); VEC row(nbCols()); putRow(iRow,row); return row;}

  /// Gets a  Vec with the copy of the column of the RealMatrix.
  /// @param iCol  number of the column
  /// @return the specified column of the matrix.
  template <typename VEC>
   VEC getColumn(std::size_t iCol) const
  {assert(iCol < nbCols()); VEC col(nbRows()); putColumn(iCol,col);return col;}

private:

  /// array data
  T m_data[N*M];
};

//////////////////////////////////////////////////////////////////////////////

/// Definition of a class RealMatrix that implements an expression template technique
/// with dynamical size
/// @author Andrea Lani
template <typename T>
class RealMatrix<T,0,0> : public ETMAT(MatExprT,RealMatrix,T,0,0) {
public:

  /// Default constructor
   RealMatrix() : ETMAT(MatExprT,RealMatrix,T,0,0)(this),
    m_owner(true), m_nrows(0), m_ncols(0), m_data(NULL) {}

  /// Constructor
   RealMatrix(std::size_t ns, std::size_t ms, T value = T()) :
    ETMAT(MatExprT,RealMatrix,T,0,0)(this), m_owner(true), m_nrows(ns), m_ncols(ms), m_data(NULL)
    {allocate();this->initialize(value);}

  // Constructor from preallocated memory
   RealMatrix(std::size_t ns, std::size_t ms, T* data) :
    ETMAT(MatExprT,RealMatrix,T,0,0)(this), m_owner(false), m_nrows(ns), m_ncols(ms), m_data(data) {}

  /// Copy constructor from expression
  template <typename EXPR>
   RealMatrix(MEETYPEV(EXPR) expr) : ETMAT(MatExprT,RealMatrix,T,0,0)(this)
  {
    m_owner = true;
    m_nrows = expr.getData()->nbRows(); m_ncols = expr.getData()->nbCols();
    allocate();
    EMATLOOP(i,0,EGETSIZER(ENMAX(0,EXPR::SIZE1)),j,0,EGETSIZEC(ENMAX(0,EXPR::SIZE2)), (*this)(i COMMA j) = expr.at(i COMMA j));
  }

  /// Copy constructor from RealMatrix
   RealMatrix(const RealMatrix<T,0,0>& orig) :
    ETMAT(MatExprT,RealMatrix,T,0,0)(this), m_owner(orig.m_owner), m_nrows(orig.m_nrows), m_ncols(orig.m_ncols)
  {
    if (m_owner) {
      allocate();
      EMATLOOP(i,0,m_nrows,j,0,m_ncols, (*this)(i COMMA j) = orig(i COMMA j));
    }
    else {m_data = orig.m_data;}
  }

  /// Default destructor
   ~RealMatrix() {free();}

  /// Overloading of assignment operator(s)
#define MAT0_EASSIGN_OP(__op__) \
  template <typename EXPR>						\
     const RealMatrix<T,0,0>& operator __op__ (MEETYPEV(EXPR) expr) \
    {EMATLOOP(i,0,EGETSIZER(ENMAX(0,EXPR::SIZE1)),j,0,EGETSIZEC(ENMAX(0,EXPR::SIZE2)), (*this)(i COMMA j) __op__ expr.at(i COMMA j)); return *this;}
MAT0_EASSIGN_OP(=)
MAT0_EASSIGN_OP(+=)
MAT0_EASSIGN_OP(-=)
MAT0_EASSIGN_OP(*=)
MAT0_EASSIGN_OP(/=)
#undef MAT0_EASSIGN_OP

  /// Overloading of assignment operator(s) with constants
#define MAT0_EASSIGN_OP_CONST(__op__) \
   const RealMatrix<T,0,0>& operator __op__ (T expr) { EVECLOOP(i, 0, size(), m_data[i] __op__ expr); return *this;}
MAT0_EASSIGN_OP_CONST(=)
MAT0_EASSIGN_OP_CONST(+=)
MAT0_EASSIGN_OP_CONST(-=)
MAT0_EASSIGN_OP_CONST(*=)
MAT0_EASSIGN_OP_CONST(/=)
#undef MAT0_EASSIGN_OP_CONST

  /// Overloading of the assignment op=() operators for RealVec
#define MAT0_EASSIGN_OP_VEC(__op__)  				\
   const RealMatrix<T,0,0>& operator __op__ (const RealVec<T,0>& diag)		\
{                                                            \
  for (std::size_t i = 0; i < m_nrows; ++i) { \
    (*this)(i,i) __op__ diag[i];     \
  }                                  \
  return *this;                      \
}
  MAT0_EASSIGN_OP_VEC(=)
  MAT0_EASSIGN_OP_VEC(+=)
  MAT0_EASSIGN_OP_VEC(-=)
  MAT0_EASSIGN_OP_VEC(*=)
  MAT0_EASSIGN_OP_VEC(/=)
#undef MAT0_EASSIGN_OP_VEC

 /// Overloading of the assignment operator "="
 /// @pre the assignee is supposed to have same size and ownership
 /// as the given RealMatrix
  const RealMatrix<T,0,0>& operator= (const RealMatrix<T,0,0>& other)
 {
   cf_assert(&other != this);
   if (m_nrows != other.m_nrows || m_ncols != other.m_ncols) {
    if (m_owner) {free();}
    m_nrows = other.m_nrows;
    m_ncols = other.m_ncols;
    if (m_owner) {allocate();}
   }
   EMATLOOP(i,0,m_nrows,j,0,m_ncols, (*this)(i COMMA j) = other(i COMMA j));
   return *this;
 }

  // Constructor from preallocated memory
   void wrap(std::size_t ns, std::size_t ms, T* data)
  {m_owner = false; m_nrows = ns; m_ncols = ms; m_data = data;}

  /// return the array size
   std::size_t size() const {return m_nrows*m_ncols;}

  /// @return a vector slice with fixed size
  template <int NV, int MV>
   RealMatrixSlice<T,NV,MV> slice(std::size_t nn, std::size_t mm)
  {return RealMatrixSlice<T,NV,MV>(&m_data[nn*m_ncols+mm], m_nrows, m_ncols);}

  /// @return a vector slice
   RealMatrixSlice<T,0,0> slice(std::size_t nn, std::size_t mm, std::size_t ns, std::size_t ms)
  {return RealMatrixSlice<T,0,0>(&m_data[nn*m_ncols+mm], m_nrows, m_ncols, ns, ms);}

  /// @return the raw data
   T* ptr() {return m_data;}

  /// return the number of rows
   std::size_t nbRows() const {return m_nrows;}

  /// return the number of columns
   std::size_t nbCols() const {return m_ncols;}

  /// Overloading of the operator"()" for assignment.
   T& operator() (std::size_t i, std::size_t j) {assert(i < m_nrows); assert(j < m_ncols); return m_data[i*m_ncols+j];}

  /// Overloading of the operator"()" (doesn't allow assignment)
   T operator() (std::size_t i, std::size_t j) const {assert(i < m_nrows); assert(j < m_ncols); return m_data[i*m_ncols+j];}

  /// Overloading of the stream operator "<<" for the output.
  /// "\n" introduced at the end of every line of the matrix.
  friend std::ostream& operator<<
    LTGT (std::ostream& out, const RealMatrix<T,0,0>& A);

  /// Overloading of the stream operator ">>" for the input
  friend std::istream& operator>>
    LTGT (std::istream& in, RealMatrix<T,0,0>& A);

  /// Accessor to individual entry
  /// @param i row index
  /// @param j column index
   T at (std::size_t i, std::size_t j) const {return (*this)(i,j);}

  /// Calculate the determinant of a 2*2 matrix
   T determ2() const;

  /// Calculate the determinant of a 3*3 matrix
   T determ3() const;

  /// Calculate the determinant of a 4*4 matrix
   T determ4() const;

  /// Invert a diagonal matrix
   void invertDiag(RealMatrix<T,0,0>& result) const;

  /// resize the array
   void resize(std::size_t n, std::size_t m, T init = T())
  {free(); m_nrows = n; m_ncols = m; allocate(); this->initialize(init);}

  /// Returns the transposes of the object matrix.
  /// The object remains untransposed.
  /// @param result transposed matrix
   void transpose(RealMatrix<T,0>& result) const
  {EMATLOOP(j, 0, m_ncols, i, 0, m_nrows, result(i COMMA j) = m_data[j*m_ncols + i]);}

  /// Puts the copy of the row of the RealMatrix in a supplied vector
  /// @param iRow  number of the row
  template <typename VEC>
   void putRow(std::size_t iRow, VEC& v) const
  {
    assert(iRow < nbRows()); assert(v.size() == nbCols());
    EVECLOOP(i, 0, nbCols(), v[i] = m_data[iRow*nbCols() + i]);
  }

  /// Puts the copy of the column of the RealMatrix in a supplied vector
  /// @param iCol  number of the column
  template <typename VEC>
   void putColumn(std::size_t iCol, VEC& v) const
  {
    assert(iCol < nbCols()); assert(v.size() == nbRows());
    EVECLOOP(i, 0, nbRows(), v[i] = m_data[nbCols()*i + iCol]);
  }

  /// Set a row of the matrix.
  /// @param row   vector with the row to set
  /// @param iRow  number of the row
  template <typename VEC>
   void setRow(const VEC& row, std::size_t iRow)
  {
    assert(row.size() == nbCols());
    std::size_t istart = iRow*nbCols();
    for (std::size_t i = 0; i < nbCols(); ++i, ++istart) {
      m_data[istart] = row[i];
    }
  }

  /// Set a column of the matrix.
  /// @param col   vector with the column to set
  /// @param iCol  number of the column
  template <typename VEC>
   void setColumn(const VEC& col, std::size_t iCol)
  {
    assert(col.size() == nbRows());
    for (std::size_t i = 0; i < nbRows(); ++i) {
      m_data[iCol+nbCols()*i] = col[i];
    }
  }

  /// Add a row of the matrix.
  /// @param row   vector with the row to set
  /// @param iRow  number of the row
  template <typename VEC>
   void addRow(const VEC& row, std::size_t iRow)
  {
    assert(row.size() == nbCols());
    std::size_t istart = iRow*nbCols();
    for (std::size_t i = 0; i < nbCols(); ++i, ++istart) {
      m_data[istart] += row[i];
    }
  }

  /// Add a column of the matrix.
  /// @param col   vector with the column to set
  /// @param iCol  number of the column
  template <typename VEC>
   void addColumn(const VEC& col, std::size_t iCol)
  {
    assert(col.size() == nbRows());
    for (std::size_t i = 0; i < nbRows(); ++i) {
      m_data[iCol+nbCols()*i] += col[i];
    }
  }

  /// Gets a Vec with the copy of the row of the RealMatrix.
  /// @param iRow  number of the row
  /// @return the specified row of the matrix.
  template <typename VEC>
   VEC getRow(std::size_t iRow) const
  {assert(iRow < nbRows()); VEC row(nbCols()); putRow(iRow,row); return row;}

  /// Gets a  Vec with the copy of the column of the RealMatrix.
  /// @param iCol  number of the column
  /// @return the specified column of the matrix.
  template <typename VEC>
   VEC getColumn(std::size_t iCol) const
  {assert(iCol < nbCols()); VEC col(nbRows()); putColumn(iCol,col);return col;}

private: // helper functions

  /// allocate the memory
   void allocate() {if (size() > 0) {assert(m_owner); m_data = new T[size()];}}

  /// free the memory
   void free() {if (m_owner && size() > 0) {delete [] m_data; m_data = NULL; m_nrows = 0; m_ncols = 0;}}

private:

  /// boolean flag to know if RealMatrix has ownership on the data
  bool m_owner;

  /// number of matrix rows
  std::size_t m_nrows;

  /// number of matrix columns
  std::size_t m_ncols;

  /// array data
  T* m_data;
};

//////////////////////////////////////////////////////////////////////////////

template <typename T, int N, int M>
std::ostream& operator<< (std::ostream& out, const RealMatrix<T,N,M>& A)
{
  for (std::size_t i = 0; i < A.nbRows(); ++i) {
    for (std::size_t j = 0; j < A.nbCols(); ++j) {
      out << A(i,j) << " " ;
    }
    out << "\n";
  }
  return out;
}

//////////////////////////////////////////////////////////////////////////////

template <typename T, int N, int M>
std::istream& operator>> (std::istream& in, RealMatrix<T,N,M>& A)
{
  for (std::size_t i = 0; i < A.nbRows(); ++i) {
    for (std::size_t j = 0; j < A.nbCols(); ++j) {
      in >> A(i,j);
    }
  }
  return in;
}

//////////////////////////////////////////////////////////////////////////////

#define MAT_DETERM2(__a__, __b__, __c__) \
template __a__ \
inline T RealMatrix<T,__b__,__c__>::determ2() const \
{						\
  assert(nbRows() == 2); assert(nbCols() == 2);		\
  return m_data[0]*m_data[3] - m_data[1]*m_data[2];	\
}
  MAT_DETERM2(<typename T COMMA int N COMMA int M>, N, M)
  MAT_DETERM2(<typename T>, 0, 0)
#undef MAT_DETERM2

//////////////////////////////////////////////////////////////////////////////

#define MAT_DETERM3(__a__, __b__, __c__) \
template __a__ \
inline T RealMatrix<T,__b__,__c__>::determ3() const \
{						\
  assert(nbRows() == 3); assert(nbCols() == 3);				\
  return m_data[0]*(m_data[4]*m_data[8] - m_data[5]*m_data[7]) -	\
    m_data[1]*(m_data[3]*m_data[8] - m_data[5]*m_data[6]) +		\
    m_data[2]*(m_data[3]*m_data[7] - m_data[4]*m_data[6]);		\
}
  MAT_DETERM3(<typename T COMMA int N COMMA int M>, N, M)
  MAT_DETERM3(<typename T>, 0, 0)
#undef MAT_DETERM3

//////////////////////////////////////////////////////////////////////////////

#define MAT_DETERM4(__a__, __b__, __c__) \
template __a__ \
T RealMatrix<T,__b__,__c__>::determ4() const \
{			      \
  assert(nbRows() == 4); assert(nbCols() == 4); \
  T d00 = m_data[0]; T d01 = m_data[1]; T d02 = m_data[2]; T d03 = m_data[3]; \
  T d04 = m_data[4]; T d05 = m_data[5]; T d06 = m_data[6]; T d07 = m_data[7]; \
  T d10d15_d14d11 = m_data[10]*m_data[15] - m_data[14]*m_data[11]; \
  T d09d15_d13d11 = m_data[ 9]*m_data[15] - m_data[13]*m_data[11]; \
  T d09d14_d13d10 = m_data[ 9]*m_data[14] - m_data[13]*m_data[10]; \
  T d08d15_d12d11 = m_data[ 8]*m_data[15] - m_data[12]*m_data[11]; \
  T d08d13_d12d09 = m_data[ 8]*m_data[13] - m_data[12]*m_data[ 9]; \
  T d08d14_d12d10 = m_data[ 8]*m_data[14] - m_data[12]*m_data[10]; \
  return  d00*(d05*(d10d15_d14d11) - d06*(d09d15_d13d11) + d07*(d09d14_d13d10)) \
    - d01*(d04*(d10d15_d14d11) - d06*(d08d15_d12d11) + d07*(d08d14_d12d10)) \
    + d02*(d04*(d09d15_d13d11) - d05*(d08d15_d12d11) + d07*(d08d13_d12d09)) \
    - d03*(d04*(d09d14_d13d10) - d05*(d08d14_d12d10) + d06*(d08d13_d12d09)); \
}
  MAT_DETERM4(<typename T COMMA int N COMMA int M>, N, M)
  MAT_DETERM4(<typename T>, 0, 0)
#undef MAT_DETERM4

//////////////////////////////////////////////////////////////////////////////

#define MAT_INVERT_DIAG(__a__, __b__, __c__) \
template __a__						\
inline void RealMatrix<T,__b__,__c__>::invertDiag(RealMatrix<T,__b__,__c__>& result) const \
{						\
  T temp = T();					\
  for (std::size_t i = 0; i < nbCols(); ++i) {	\
    temp = (*this)(i,i);			\
    assert(std::abs(temp) > 0.);		\
    result(i,i) = 1./temp;			\
  }						\
}
  MAT_INVERT_DIAG(<typename T COMMA int N COMMA int M>, N, M)
  MAT_INVERT_DIAG(<typename T>, 0, 0)
#undef MAT_INVERT_DIAG

//////////////////////////////////////////////////////////////////////////////

}   /*--- End of namespace Common ---*/

//////////////////////////////////////////////////////////////////////////////

#endif
