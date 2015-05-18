// =========================================================================
/* Copyright 2011 Sebastien Bougleux
   
   This file is part of GeoComp.
   
   GeoComp is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.
   
   GeoComp is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU Lesser General Public License for more details.
   
   You should have received a copy of the GNU General Public License,
   and a copy of the GNU Lesser General Public License, along with 
   GeoComp. If not, see <http://www.gnu.org/licenses/>.
*/
// =========================================================================
/**
 * @file polynomials.hh
 */
// =========================================================================
#ifndef __POLYNOMIALS_HH__
#define __POLYNOMIALS_HH__

#include <cmath>
#include <limits>

/**
 * @namespace GeoComp
 * @brief Tools for computational and numerical geometry.
 */
namespace GeoComp
{
  // ----------------------------------------------------------------
  template <typename FT>
  struct Complex
  {
    FT re;
    FT img;
  };
  // ----------------------------------------------------------------
  /**
   * @brief Compute the real solution(s) of the quadratic equation at^2+bt+c=0.
   * @param[in] a value of a.
   * @param[in] b value of b.
   * @param[in] c value of c.
   * @param[out] t1 1st solution, involving -sqrt(delta).
   * @param[out] t1 1st solution, involving sqrt(delta).
   * @return the discriminant delta.
   */
  template <typename FT>
  double quadratic(FT a, FT b, FT c, FT &t1, FT &t2)
  {
    FT delta = b*b - 4.0*a*c;
    if (delta > 0.0)
    {
      t1 = (-b - std::sqrt(delta)) / (2.0*a);
      t2 = (-b + std::sqrt(delta)) / (2.0*a);
    }
    else
    {
      if (delta < 0.0) return delta;
      t2 = 0;
      t1 = -b / (2.0*a);
    }
    return delta;
  }
  // ----------------------------------------------------------------
  /**
   * @brief Compute the complex solution(s) of the quadratic equation at^2+bt+c=0.
   * @param[in] a value of a.
   * @param[in] b value of b.
   * @param[in] c value of c.
   * @param[out] z1 1st solution, involving -sqrt(delta).
   * @param[out] z1 1st solution, involving sqrt(delta).
   * @return the discriminant delta.
   */
  template <typename FT>
  double quadratic(FT a, FT b, FT c, Complex<FT> &z1, Complex<FT> &z2)
  {
    FT delta = b*b - 4.0*a*c;
    if (delta > 0.0)
    {
      z1.img = z2.img = 0.0;
      z1.re = (-b - std::sqrt(delta)) / (2.0*a);
      z2.re = (-b + std::sqrt(delta)) / (2.0*a);
    }
    else
    {
      if (delta < 0.0)
      {
	z1.re = z2.re = -b / (2.0*a);
	z1.img = -std::sqrt(-delta) / (2.0*a);
	z2.img = -z1.img;
      }
      else // delta == 0.0
      {
	z1.img = z2.img = 0;
	z1.re = z2.re = -b / (2.0*a);
      }
    }
    return delta;
  }
  // ----------------------------------------------------------------
  template <typename FT>
  FT cubicRoot(FT x)
  {
    return (x < 0.0 ?
	    -std::exp((1.0/3.0)*std::log(-x)) :
	    std::exp((1.0/3.0)*std::log(x)));
  }
  // ----------------------------------------------------------------
  template <typename FT>
  FT cubicRoot(FT re, FT img)
  {
    FT mod, arg, re_m, img_m, z_re;
    mod = sqrt(re * re + img * img);
    re_m = re / mod;
    img_m = img / mod;
    if (re_m < 0.0 && img_m > 0.0) arg = std::acos(re_m);
    else arg = std::asin(img_m);
    z_re = std::exp((1.0/3.0)*std::log(mod)) * std::cos(arg/3.0);
    return z_re;
  }
  // ----------------------------------------------------------------
  template <typename FT>
  FT cubic(FT a, FT b, FT c, FT d, Complex<FT> &z1, Complex<FT> &z2, Complex<FT> &z3)
  {
    FT delta, p, q, s;
    b /= a; c /= a; d /= a; a = 1.0;
    s = -b/3.0;
    p = -b*b/3.0 + c;
    q = (b/27.0) * (2.0*(b*b) - 9.0*c) + d;
    delta = q*q + (4.0/27.0) * (p*p*p);
    if (delta > 0.0)
    {
      z1.img = 0.0;
      z1.re = s + cubicRoot((-q+std::sqrt(delta))/2.0) + cubicRoot((-q-std::sqrt(delta))/2.0);
      quadratic(1.0, b+z1.re, c+z1.re*(b+z1.re), z2, z3);
    }
    else
    {
      if (delta < 0.0)
      {
	z1.img = z2.img = z3.img = 0.0;
	z1.re = cubicRoot(-q/2.0, std::sqrt(-delta)/2.0);
	z1.re += z1.re + s;
	quadratic(1.0, b+z1.re, c+z1.re*(b+z1.re), z2, z3);
      }
      else
      {
	z1.img = z2.img = z3.img = 0.0;
	z1.re = s + 3.0*(q/p);
	z2.re = z3.re = s + -(3.0/2.0)*(q*p);
      }
    }
    return delta;
  }
  // ----------------------------------------------------------------
  template <typename FT>
  bool quartic(FT a, FT b, FT c, FT d, FT e, Complex<FT> &z0, Complex<FT> &z1, Complex<FT> &z2, Complex<FT> &z3)
  {
    Complex<FT> a0, b0;
    FT s, p, q, r, x0;
    e /= a;
    d /= a;
    c /= a;
    b /= a;
    a = 1.0;
    s = -b / 4.0;
    p = (-3.0 * b * b) / 8.0 + c;
    q = (b / 2.0 * b / 2.0 * b / 2.0) - 0.5 * b * c + d;
    r = -3.0 * (b / 4.0 * b / 4.0 * b / 4.0 * b / 4.0) + c * (b / 4.0 * b / 4.0) - 0.25 * b * d + e;
    cubic(8.0, -(4.0 * p), -(8.0 * r), 4.0 * r * p - (q * q), z0, z1, z2);
    x0 = z0.re;
    a0.re = -p + 2.0 * x0;
    if (a0.re < 0.0)
    {
      // TODO
      return false;
    } 
    else
    {
      a0.img = b0.img = 0.0;
      a0.re = sqrt(a0.re);
      b0.re = -q / (2.0 * a0.re);
      quadratic(1.0, -a0.re, x0 - b0.re, z0, z1);
      quadratic(1.0, a0.re, x0 + b0.re, z2, z3);
      z0.re += s;
      z1.re += s;
      z2.re += s;
      z3.re += s;
    }
    return true;
  }
  // ----------------------------------------------------------------
  // return at^3+bt^2+ct+d
  template <typename FT>
  inline FT cubicPoly(const FT &t, const FT &a, const FT &b, const FT &c, const FT &d)
  { return a*(t*t*t) + b*(t*t) + c*t + d; }
  // ----------------------------------------------------------------
  // return at^4+bt^3+ct^2+dt+e
  template <typename FT>
  inline FT quarticPoly(const FT &t, const FT &a, const FT &b, const FT &c, const FT &d, const FT &e)
  { FT t2 = t*t; return a*(t2*t2) + b*(t2*t) + c*t2 + d*t + e; }
  // ----------------------------------------------------------------
  template <typename FT>
  FT quartic_newton(const FT &a, const FT &b, const FT &c, const FT d, const FT &e, 
		    FT sol, FT epsilon = 1.0e-10, int nbItMax = 10)
  {
    FT ad = 4.0*a, bd = 3.0*b, cd = 2.0*c, error = std::numeric_limits<FT>::max(), tmp;
    for (; nbItMax > 0 && error > epsilon; nbItMax--)
    {
      tmp = sol - quarticPoly(sol,a,b,c,d,e) / cubicPoly(sol,ad,bd,cd,d);
      error = std::abs(sol-tmp);
      sol = tmp;
    }
    return sol;
  }

} // end namespace

#endif
