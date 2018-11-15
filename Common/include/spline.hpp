#ifndef SU2_SPLINE_HPP
#define SU2_SPLINE_HPP

#include "datatype_structure.hpp"
#include "datatypes/vectorT.hpp"

#include <vector>

namespace MathTools {
  
  using RealVec = Common::RealVec;

  /*!
    * \brief Sets the second derivative coefficients for natural spline.
    * \param[in] x The vector of values to be tabulated (input)
    * \param[in] y The vector of values tabulated,i.e yi = f(xi) (input)
    * \param[in] yp1 Value of first derivative in the first point (input)
    * \param[in] ypn Value of first derivative in the last point (input)
    * \param[out] y2 Coefficients of second derivative for interpolating function (output)
  */
  void SetSpline(RealVec& x, const RealVec& y, const su2double yp1, const su2double ypn, RealVec& y2);

  /*!
    * \brief Sets the second derivative coefficients for natural spline.
    * \param[in] x The vector of values to be tabulated (input)
    * \param[in] y The vector of values tabulated,i.e yi = f(xi) (input)
    * \param[in] y2 Coefficients of second derivative for interpolating function (output)
    * \param[in] x The value we want to tabulate through spline (input)
  */
  su2double GetSpline(RealVec& x, const RealVec& y, const RealVec& y2, su2double value);
}

#endif
