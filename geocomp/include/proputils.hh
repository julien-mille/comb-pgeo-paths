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
 * @file proputils.hh
 */
// =========================================================================
#ifndef __PROPUTILS_HH__
#define __PROPUTILS_HH__

#include <vector>
#include <algorithm>

/**
 * @namespace GeoComp
 * @brief Tools for computational and numerical geometry.
 */
namespace GeoComp
{
  // ----------------------------------------------------------------
  /** @class SaddlePair
   */
  // ----------------------------------------------------------------
  template <typename IT, typename FT>
  struct SaddlePair
  {
    // ----------------------------------------------------------------
    IT lab1;  // label of the 1st node
    IT lab2;  // label of the 2nd node
    // ----------------------------------------------------------------
    IT idx1;  // index of the saddle associated to lab1
    IT idx2;  // index of the saddle associated to lab2
    // ----------------------------------------------------------------
    FT d;
    // ----------------------------------------------------------------
    SaddlePair(IT label1, IT index1, IT label2, IT index2, FT dst)
      : lab1(label1), lab2(label2), idx1(index1), idx2(index2), d(dst) { }
    // ----------------------------------------------------------------
    ~SaddlePair() { }
    // ----------------------------------------------------------------
    bool update(const IT &label1, const IT &index1, const IT &label2, const IT &index2, const FT &dst)
    {
      if (dst < d)
      {
	lab1 = label1; lab2 = label2; idx1 = index1; idx2 = index2; d = dst;
	return true;
      }
      return false;
    }
    // ----------------------------------------------------------------
  };
  // ----------------------------------------------------------------
  /** @class SaddlePairCompare
   */
  // ----------------------------------------------------------------
  class SaddlePairCompare
  {
  public:
    // ----------------------------------------------------------------
    SaddlePairCompare() { }
    // ----------------------------------------------------------------
    ~SaddlePairCompare() { }
    // ----------------------------------------------------------------
    template <typename IT, typename FT>
    bool operator()(const SaddlePair<IT,FT> &s1, const SaddlePair<IT,FT> &s2)
    { return (s1.d < s2.d); }
    // ----------------------------------------------------------------
  };
  // ----------------------------------------------------------------
  /** @class VoronoiEdge
   */
  // ----------------------------------------------------------------
  template <typename IT, typename FT>
  struct VoronoiEdge
  {
    typedef std::vector<SaddlePair<IT,FT> > VPairs;
    // ----------------------------------------------------------------
    std::vector<SaddlePair<IT,FT> > vpairs;
    // ----------------------------------------------------------------
    VoronoiEdge(IT label1, IT index1, IT label2, IT index2, FT dst)
    { vpairs.push_back(SaddlePair<IT,FT>(label1,index1,label2,index2,dst)); }
    // ----------------------------------------------------------------
    ~VoronoiEdge() { }
    // ----------------------------------------------------------------
    bool update(const IT &label1, const IT &index1, const IT &label2, const IT &index2, const FT &dst)
    {
      vpairs.push_back(SaddlePair<IT,FT>(label1,index1,label2,index2,dst));
      return true;
    }
    // ----------------------------------------------------------------
    void sort()
    {
      SaddlePairCompare spc;
      std::sort(vpairs.begin(),vpairs.end(),spc);
    }
    // ----------------------------------------------------------------
  };
  // ----------------------------------------------------------------
  /** @class DelGraphArray
   *  @brief Array-based Delaunay graph (edges encoded by a node-node static matrix)
   */
  // ----------------------------------------------------------------
  template <typename IT, typename FT, class Edg = SaddlePair<IT,FT> >
  class DelGraphArray
  {
  public:
    typedef Edg Edge;
  protected:
    IT _nb_nodes;
    Edge ***_g;
  public:
    // ----------------------------------------------------------------
    DelGraphArray() : _nb_nodes(0), _g(NULL) { }
    // ----------------------------------------------------------------
    DelGraphArray(IT nb_nodes) : _nb_nodes(0), _g(NULL)
    { reserve(nb_nodes); }
    // ----------------------------------------------------------------
    ~DelGraphArray()
    { clear(); }
    // ----------------------------------------------------------------
    void clear()
    {
      if (_g)
      {
	for (IT i = 0, j = 0; i < _nb_nodes; i++)
	{
	  for (j = i+1; j < _nb_nodes; j++)
	    if (_g[i][j]) { delete _g[i][j]; _g[j][i] = _g[i][j] = NULL; }
	  delete[] _g[i];
	}
	delete[] _g;
	_g = NULL;
	_nb_nodes = 0;
      }
    }
    // ----------------------------------------------------------------
    void reserve(IT nb_nodes)
    {
      if (nb_nodes != _nb_nodes)
      {
	clear();
	_nb_nodes = nb_nodes;
	_g = new Edge**[nb_nodes];
	for (IT i = 0, j = 0; i< nb_nodes; i++)
	{
	  _g[i] = new Edge*[nb_nodes];
	  for (j = 0; j < nb_nodes; j++) _g[i][j] = NULL;
	}
      }
      else
      {
	for (IT i = 0, j = 0; i< _nb_nodes; i++)
	  for (j = i+1; j < nb_nodes; j++)
	  {
	    if (_g[i][j])
	    { 
	      delete _g[i][j];
	      _g[j][i] = _g[i][j] = NULL;
	    }
	  }
      }
    }
    // ----------------------------------------------------------------
    template <class Container>
    void edges(Container &le)
    {
      for (IT i = 0, j = 0; i < _nb_nodes; i++)
	for (j = i+1; j < _nb_nodes; j++)
	  if (_g[i][j]) le.push_back(_g[i][j]);
    }
    // ----------------------------------------------------------------
    bool newEdge(const IT &lab1, const IT &lab2, const IT &idx1, const IT &idx2, FT d)
    {
      Edge *e = _g[lab1][lab2];
      if (e == NULL)
      {
	_g[lab2][lab1] = _g[lab1][lab2] = new Edge(lab1,idx1,lab2,idx2,d);
	return true;
      }
      else e->update(lab1,idx1,lab2,idx2,d);
      return false;
    }
    // ----------------------------------------------------------------
  
  };
  
} // end namespace

#endif
