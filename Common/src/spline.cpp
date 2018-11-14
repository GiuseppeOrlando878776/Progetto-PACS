#include "../include/spline.hpp"

#include <algorithm>
#include <exception>
#include "../include/su2_assert.hpp"

namespace MathTools {
  //
  /*--- Sets the second derivative coefficients for natural spline ---*/
  void SetSpline(RealVec& x, const RealVec& y, const su2double yp1, const su2double ypn, RealVec& y2) {

    auto n = x.size();
    SU2_Assert(n > 1,"You have only one datum: cannot generate spline");
    SU2_Assert(y.size() == n,"The dimension of vectors x and y from data are not the same");
    SU2_Assert(std::is_sorted(x.cbegin(),x.cend()),"The vector of x is not ordered");
    auto it = std::unique(x.begin(), x.end());
    SU2_Assert(it == x.end(),"Some values that have been tabulated are not unique");

    unsigned long i, k;
    double p, qn, sig, un;
    RealVec u(n);
    y2.clear();
    y2.reserve(n);

    if (yp1 > 0.99e30) 			// The lower boundary condition is set either to be "natural"
      y2.push_back(0.0);
    else {									// or else to have a specified first derivative.
      y2.push_back(-0.5);
      u[0]=(3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
    }

    for (i=2; i<n; ++i) {									//  This is the decomposition loop of the tridiagonal algortihm
      sig = (x[i-1]-x[i-2])/(x[i]-x[i-2]);		//	y2 and u are used for temporary
      p = sig*y2[i-2] + 2.0;										//	storage of the decomposed
      y2.push_back((sig-1.0)/p);									//	factors.
      u[i-1] = (y[i]-y[i-1])/(x[i]-x[i-1]) - (y[i-1]-y[i-2])/(x[i-1]-x[i-2]);
      u[i-1] = (6.0*u[i-1]/(x[i]-x[i-2])-sig*u[i-2])/p;
    }

    if (ypn > 0.99e30)						// The upper boundary condition is set either to be
      qn=un=0.0;									// "natural"
    else {												// or else to have a specified first derivative.
      qn=0.5;
      un=(3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
    }
    y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0);
    for (k=n-1; k>0; --k)					// This is the backsubstitution loop of the tridiagonal
      y2[k-1]=y2[k-1]*y2[k]+u[k-1];	  // algorithm.

  }

  //
  /*--- Gets the value from natural spline ---*/
  su2double GetSpline(RealVec& x, const RealVec& y, const RealVec& y2, su2double value) {

    auto n = x.size();
    SU2_Assert(n > 1,"You have only one datum: cannot use spline");
    SU2_Assert(y.size() == n,"The dimension of vectors x and y from data are not the same");
    SU2_Assert(y2.size() == n,"The dimension of vector of coefficient of second derivative doesn't match");
    SU2_Assert(std::is_sorted(x.cbegin(),x.cend()),"The vector of x is not ordered");
    auto it = std::unique(x.begin(), x.end());
    SU2_Assert(it == x.end(),"Some values that have been tabulated are not unique");

    bool is_equispaced = true;
    const auto step = x[1] - x[0];
    for(std::size_t i=2; i<n and is_equispaced; ++i) {
      const auto curr_step = x[i] - x[i-1];
      if(curr_step != step)
        is_equispaced = false;
    }
    SU2_Assert(is_equispaced,"This version of spline works only with equispaced data since it uses integer division");

    if(value < x[0] or value > x[n-1])
      throw std::out_of_range("The required temperature is out of data range");

    unsigned long klo, khi;
    double h, b, a, result;
    klo = (value - x[0])/step + 1;
    khi = klo + 1;
    h = step;
    a = (x[khi-1] - value)/h;
    b = (value - x[klo-1])/h;
    // Cubic spline polynomial is now evaluated.
    result = a*y[klo-1] + b*y[khi-1] + ((a*a*a-a)*y2[klo-1] + (b*b*b-b)*y2[khi-1])*(h*h)/6.0;

    return result;

  }

  }
