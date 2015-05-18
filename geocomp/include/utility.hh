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
 * @file utility.hh
 */
// =========================================================================
#ifndef __UTILITY_HH__
#define __UTILITY_HH__

#include <cstddef>
#include <fstream>
#include <vector>
#include <cmath>

/**
 * @namespace GeoComp
 * @brief Tools for computational and numerical geometry.
 */
namespace GeoComp
{
  // ----------------------------------------------------------------
  /**
   * @brief Test if an integer is odd.
   * @param[in] k Value of an integer.
   * @return true if k is odd, false else.
   */
  #define ODD(k) (k & 1)
  // ----------------------------------------------------------------
  template <typename T>
  inline T scalarProduct(const T &v1x, const T &v1y, const T &v2x, const T &v2y)
  {
    return (v1x*v2x + v1y*v2y);
  }
  // ----------------------------------------------------------------
  template <typename T>
  inline T norm(const T &vx, const T &vy) { return std::sqrt(vx*vx+vy*vy); }
  // ----------------------------------------------------------------
  template <typename T>
  inline T squaredNorm(const T &vx, const T &vy) { return (vx*vx+vy*vy); }
  // ----------------------------------------------------------------
  template <typename T>
  inline void normalize(T &vx, T &vy)
  {
    T nv = norm(vx,vy);
    vx /= nv;
    vy /= nv;
  }
  // ----------------------------------------------------------------
  template <typename FT>
  void path2svg(std::vector<std::pair<FT,FT> > &path, const char *filename)
  {
    std::ofstream ofile(filename);
    if (!ofile) return;
    ofile << "<!DOCTYPE html>\n" << "<html>\n" << "<body>\n"
	  << "<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\">\n";
    //std::vector<std::pair<FT,FT> >::iterator eit = path.begin(), eend = path.end();
    for (long i = 0; i < path.size()-1; i++)
      ofile << "<line x1=\"" << path[i].second << "\" y1=\"" << path[i]->first
	    << "\" x2=\"" << path[i+1]->second << "\" y2=\""
	    << path[i+1]->first << "\" style=\"stroke:rgb(255,0,0);stroke-width:2\"/>\n";
    ofile << "</svg>\n";
    ofile << "</body>\n</html>\n";
    ofile.close();
  }
} // end namespace

#endif
