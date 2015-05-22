// =========================================================================
/* Copyright 2012 Sebastien Bougleux

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
 * @file grids.hh
 */
// =========================================================================
#ifndef __GRIDS_HH__
#define __GRIDS_HH__

#include <limits>
#include <iostream>
#include <cmath>
#include <vector>
#include <queue>
#include <cstring>
#include "utility.hh"

namespace GeoComp
{
  /**
   * @brief Loop to visit all interior nodes.
   * @param[in] g 2D grid (RectGrid).
   * @param[in] i Index of the node.
   * @param[in] IT Type of the index (must be integer).
   * @par Example
   * @code
   RectGrid<long> g(250,250);
   double *fct = new double[g.size()];
   long i;
   for_grid_interior(g,i,long)
      fct[i] = someFunction();
   *  @endcode
   */
#define for_grid_interior(g,i,IT)					\
  for (IT n = g.firstInteriorBlockBegin(), stop = g.firstInteriorBlockEnd(); \
       stop < g.lastInteriorBlockEnd(); stop = g.nextInteriorBlockEnd(stop), n = g.nextInteriorBlockBegin(n)) \
    for (i = n; i < stop; i++)
  // ------------------------------------------------------------
#define for_interior_nodes(i,IT)					\
  for (IT n = firstInteriorBlockBegin(), stop = firstInteriorBlockEnd(), i; \
       stop < lastInteriorBlockEnd(); stop = nextInteriorBlockEnd(stop), n = nextInteriorBlockBegin(n)) \
    for (i = n; i < stop; i++)
  // ------------------------------------------------------------
#define FOREACH_4_NEIGHBOR(k) for(k=0;k<4;k++)
#define FOREACH_8_NEIGHBOR(k) for(k=0;k<7;k++)
#define FOREACH_4_8_NEIGHBOR(k) for(k=0;k<7;k+=2)
#define FOREACH_8DIAG_NEIGHBOR(k) for(k=1;k<8;k+=2)
  // ------------------------------------------------------------
  /**
   * @class Grid
   * @brief 2D rectangular grid.
   */
  template <typename IT = long>
  class Grid
  {
  protected:
    // ------------------------------------------------------------
    IT _size;
    IT _rows;
    IT _cols;
    // ------------------------------------------------------------
    IT _n4[4];
    IT _n8[8];
    IT *_n[2];
    static const unsigned short _opposite4[4];
    static const unsigned short _opposite8[8];
    static const unsigned short *_opposite[2];
    // ------------------------------------------------------------
    void constructNeighbors()
    {
      _n4[0] = _n8[0] = 1;
      _n4[1] = _n8[2] = _cols;
      _n4[2] = _n8[4] = -1;
      _n4[3] = _n8[6] = -_cols;
      _n8[1] = _cols+1;
      _n8[3] = _cols-1;
      _n8[5] = -_cols-1;
      _n8[7] = -_cols+1;
      _n[0] = _n4; _n[1] = _n8;
    }
    // ------------------------------------------------------------
    inline unsigned short nodeValence(unsigned short ringRadius)
    { return (ringRadius == 0 ? 4 : 8); }
    // ------------------------------------------------------------
  public:
    typedef IT IndexType;
    // ------------------------------------------------------------
    /**
     * @brief Constructor from grid dimensions.
     * @param[in] r Number of rows of the grid.
     * @param[in] c Number of columns of the grid.
     */
    Grid(IT r, IT c) : _size(r*c), _rows(r), _cols(c) { _n[0] = _n[1] = NULL; constructNeighbors(); }
    // ------------------------------------------------------------
    /**
     * @brief Destructor.
     */
    ~Grid() { }
    // ------------------------------------------------------------
    /**
     * @brief Size of the grid (number of nodes).
     */
    inline const IT& size() const { return _size; }
    // ------------------------------------------------------------
    /**
     * @brief Number of nodes on a row of the grid.
     */
    inline const IT& rows() const { return _rows; }
    // ------------------------------------------------------------
    /**
     * @brief Number of nodes on a column of the grid.
     */
    inline const IT& cols() const { return _cols; }
    // ------------------------------------------------------------
    /**
     * @brief Index of a node (its offset), provided its position (row,colunm).
     * @param[in] r Row of the node.
     * @param[in] c Column of the node.
     * @note No validity check is performed on the position values (r,c).
     */
    inline IT index(IT r, IT c) const { return c + r * _cols; }
    // ------------------------------------------------------------
    /*inline IT index(IT *r, IT *c, IT nb_elements, IT *idx)
    {
      for (IT k = 0; k < nb_elements; k++)
	idx[k] = (_rows_or_cols ? c[k] + r[k] * _subsize : r[k] + c[k] * _subsize);
	}*/
    // ------------------------------------------------------------
    /**
     * @brief Position (row,column) of a node provided its index.
     * @param[in] idx Index of the node (offset).
     * @param[out] r Corresponding row.
     * @param[out] c Corresponding columns.
     * @note No validity check is performed on the value of idx.
     */
    inline void point(const IT &idx, IT &r, IT &c) const
    { r = idx / _cols; c = idx - r*_cols; }
    // -------------------------------------------------------------------
    template <class Container>
    inline void points(const Container &p_idx, IT *p_rc) const
    {
      IT r, c;
      for (IT i = 0, j = 0; i < p_idx.size(); i++, j+=2)
      {
	point(p_idx[i],r,c);
	p_rc[j] = r;
	p_rc[j+1] = c;
      }
    }
    // -------------------------------------------------------------------
    /**
       @brief Row associated to a given floatting-point value
       @param[in] y Floatting-point y-coordinate.
       @return -1 if (y<-0.5 or y>rows()-0.5), 0 if y=0.5, round(y) else.
     */
    template <typename FT>
    IT row(FT y) const
    {
      if (y < -0.5) return -1;
      if (y == -0.5) return 0;
      if (y > (FT)(_rows)-(FT)0.5) return -1;
      if (y == (FT)(_rows-1)+(FT)0.5) return (_rows-1);
      return (IT)(y + (FT)0.5);
    }
    // -------------------------------------------------------------------
    /**
       @brief Column associated to a given floatting-point value
       @param[in] x Floatting-point x-coordinate.
       @return -1 if (y<-0.5 or y>rows()-0.5), 0 if y=0.5, round(y) else.
     */
    template <typename FT>
    IT col(FT x) const
    {
      if (x < -0.5) return -1;
      if (x == -0.5) return 0;
      if (x > (FT)(_cols)-0.5) return -1;
      if (x == (FT)(_cols-1)+0.5) return (_cols-1);
      return (IT)(x + (FT)0.5);
    }
    // -------------------------------------------------------------------
    template <typename ET>
    void row2col(ET *f) const
    {
      IT i = 0, r, c;
      ET *ft = new ET[_size];
      copy(f,ft);
      for (; i < _size; i++)
      {
	c = i / _rows;
	r = i - c*_rows;
	f[i] = ft[index(r,c)];
      }
      delete[] ft;
    }
   // -------------------------------------------------------------------
    template <typename ET>
    void col2row(ET *f) const
    {
      IT i = 0, r, c;
      ET *ft = new ET[_size];
      copy(f,ft);
      for (; i < _size; i++)
      {
	r = i / _cols;
	c = i - r*_cols;
	f[i] = ft[r + c * _rows];
      }
      delete[] ft;
    }
    // ------------------------------------------------------------
    template <typename FT>
    void labelize(FT *f, IT *l)
    {
      std::queue<IT> q;
      IT lab = 1, j, k, m;
      fill(l,-1);
      for (IT i = 0; i < _size; i++)
      {
	if (f[i] > 0)
	{
	  if (l[i] != -1) continue;
	  q.push(i);
	  while (!q.empty())
	  {
	    j = q.front(); q.pop();
	    l[j] = lab;
	    for (k = 0; k < 8; k++)
	    {
	      m = neighbor8(j,k);
	      if (l[m] != -1) continue;
	      q.push(m);
	    }
	  }
	  lab++;
	}
	else l[i] = 0;
      }
    }
    // ------------------------------------------------------------
    /**
     * @brief Provide the neighborhood index associated to a given adjacency.
     * @param[in] adjacency Adjacency.
     * @return Neighborhood index (0 if adjacency==4, 1 else).
     */
    inline unsigned short nodeNeighborhoodIndex(unsigned short adjacency) const
    { return (adjacency == 4 ? 0 : 1); }
    // ------------------------------------------------------------
    /**
     * @brief Index of the first interior node.
     */
    inline IT firstInteriorBlockBegin() const { return (_cols + 1); }
    // ------------------------------------------------------------
    /**
     * @brief Index of the first interior node.
     */
    inline IT firstInteriorRowFirstNode() const { return (_cols + 1); }
    // ------------------------------------------------------------
    /**
     * @brief Index of the first interior node.
     */
    inline IT firstInteriorRowLastNode() const { return (2*_cols - 2); }
    // ------------------------------------------------------------
    /**
     * @brief Index of the first interior node.
     */
    inline IT lastInteriorRowFirstNode() const { return (_size - 2*_cols + 2); }
    // ------------------------------------------------------------
    /**
     * @brief Index of the first interior node.
     */
    inline IT lastInteriorRowLastNode() const { return (_size - _cols - 2); }
    // ------------------------------------------------------------
    /**
     * @brief Index of the last interior node of the first block (row or column).
     */
    inline IT firstInteriorBlockEnd() const { return (2*_cols-1); }
    // ------------------------------------------------------------
    /**
     * @brief Index of the (last+1) interior node.
     */
    inline IT lastInteriorBlockEnd() const { return (_size-_cols); }
    // ------------------------------------------------------------
    /**
     * @brief Index of the last interior node of the next block (row or column).
     * @param[in] i Index of the last interior node of a block.
     */
    inline IT nextInteriorBlockEnd(const IT &i) const { return (i+_cols); }
    // ------------------------------------------------------------
    /**
     * @brief Index of the first interior node of the next block (row or column).
     * @param[in] i Index of the first interior node of a block.
     */
    inline IT nextInteriorBlockBegin(const IT &i) const { return (i+_cols); }
    // ------------------------------------------------------------
    /**
     * @brief Index of the k-st neighbor of a node in 4-adjacency.
     * @param[in] idx Index of a node.
     * @param[in] k Neighbor of idx, which must be in [0,3].
     * @return Index of the k-st neighbor of a node in 4-adjacency.
     */
    inline IT neighbor4(const IT idx, const unsigned short k) const { return (idx + _n4[k]); }
    // ------------------------------------------------------------
    /**
     * @brief Index of the k-st neighbor of a node in 8-adjacency.
     * @param[in] idx Index of a node.
     * @param[in] k Neighbor of idx, which must be in [0..7].
     * @return Index of the k-st neighbor of a node in 8-adjacency.
     */
    inline IT neighbor8(const IT &idx, const unsigned short &k) const { return (idx + _n8[k]); }
    // ------------------------------------------------------------
    /**
     * @brief Index of the k-st neighbor of a node in 4 or 8-adjacency.
     * @param[in] idx Index of a node.
     * @param[in] radius Adjacency, 0 => 4-adjaceny, 1 => 8-adjacency.
     * @param[in] k Neighbor of idx, which must be in [0..3] in 4-adjacency, and in [0..7] in 8-adjacency.
     * @return Index of the k-st neighbor of a node.
     */
    inline IT neighbor(IT idx, unsigned short radius, unsigned short k) const
    { return (idx + _n[radius][k]); }
    // ------------------------------------------------------------
    /**
     * @brief Opposite neighbor index in 4-adjacency.
     * @param[in] k Neighbor index of a node (in [0..3]).
     * @return Opposite neighbor index.
     */
    inline const unsigned short& opposite4(const unsigned short &k) const
    { return _opposite4[k]; }
    // ------------------------------------------------------------
    /**
     * @brief Opposite neighbor index in 8-adjacency.
     * @param[in] k Neighbor index of a node (in [0..7]).
     * @return Opposite neighbor index.
     */
    inline const unsigned short& opposite8(const unsigned short &k) const
    { return _opposite8[k]; }
    // ------------------------------------------------------------
    /**
     * @brief Opposite neighbor index in 8-adjacency.
     * @param[in] k Neighbor index of a node (in [0..7]).
     * @return Opposite neighbor index.
     */
    inline const unsigned short& opposite(const unsigned short &radius, const unsigned short &k) const
    { return _opposite[radius][k]; }
    // ------------------------------------------------------------
    /**
     * @brief Test if a node given by its index is on the grid boundary.
     * @param[in] i Index of the node on the grid.
     * @return true if the node is on the grid boundary, false else.
     */
    inline bool onBoundary(const IT &i) const
    { return ((i < _cols) || (_size < _cols + 1 + i) || (i % _cols == 0) ||
	      ((i+1) % _cols == 0)); }
    // ------------------------------------------------------------
    template <typename FT>
    bool isLocalMax(FT *fi)
    {
      unsigned short count = 0;
      for (unsigned short k = 0; k < 8; k++)
	if (fi[_n8[k]] > *fi) return false;
	else if (fi[_n8[k]] == *fi) count++;
      return (count == 8 ? false : true);
    }
    // ------------------------------------------------------------
    template <typename FT, class Container>
    void localMax(FT *f, Container &lm)
    {
      IT r,c;
      fillBoundary(f,std::numeric_limits<FT>::min());
      for_interior_nodes(i,IT) if (isLocalMax(&f[i])) lm.push_back(i);
	// {
	//   point(i,r,c);
	//   std::cerr << "[" << r << ";" << c << "],";
	// }
    }
    // ------------------------------------------------------------
    template <typename FT, typename IT2>
    void voronoiPartition(FT *D1, FT *D2, IT2 *V)
    {
      fillBoundary(V,0);
      for_interior_nodes(i,IT)
	if (D1[i] < D2[i]) V[i] = 1;
	else
	  if (D1[i] > D2[i]) V[i] = 2;
	  else V[i] = 3;
    }
    // ------------------------------------------------------------
    template <typename FT, typename LT>
    void saddles(FT *D1, FT *D1r, FT *D1c, FT *D2r, FT *D2c, FT *D2, LT *V, FT *VE)
    {
      unsigned short k;
      IT j;
      FT nrm1, nrm2, dsr, dsc;
      bool is_max, is_ve;
      voronoiPartition(D1,D2,V);
      fill(VE,0);
      for_interior_nodes(i,IT)
      {
	is_ve = false;
	is_max = true;
	for (k = 0; k < 8; k++)
	{
	  j = neighbor8(i,k);
	  if (V[i] != V[j] && V[j] > 0) is_ve = true;
	  if (V[i] == 1 && D1[j] > D1[i]) is_max = false;
	  else if (V[i] == 2 && D2[j] > D2[i]) is_max = false;
	}
	if (is_ve)
	{
	  if (is_max)
	  nrm1 = std::sqrt(D1r[i]*D1r[i] + D1c[i]*D1c[i]);
	  nrm2 = std::sqrt(D2r[i]*D2r[i] + D2c[i]*D2c[i]);
	  dsr = D1r[i] / nrm1 + D2r[i] / nrm2;
	  dsc = D1c[i] / nrm1 + D2c[i] / nrm2;
	  VE[i] = std::sqrt(dsr*dsr + dsc * dsc);
	}

      }
    }
    // ------------------------------------------------------------
    /**
     * @brief Nearest-neighbors interpolation of a function at a given floatting point coordinates, knowing the value of the function on the gird nodes (integer coordinates).
     * @param[in] pr Floatting point coordinate on rows.
     * @param[in] pc Floatting point coordinate on columns.
     * @param[in] f Function on grid nodes.
     * @param[in] outside Value representing the outside of the definition domain of f. Nodes of the grid marked with outside (f[i]=outside) are not taken into account in the interpolation process. Usually it is set to the max or the min value of the type FT.
     */
    template <typename FT>
    FT nnInterpolation(FT pr, FT pc, const FT *f, FT outside)
    {
      IT pri = (IT)pr, pci = (IT)pc;
      FT res = 0, tmp, w = 0;
      // --- special cases ------------------------
      if (pr-(FT)pri == 0) // same row
      {
	if (pc-(FT)pci == 0) return f[index(pri,pci)]; // integer point
	if (f[index(pri,pci)] != outside)
	{
	  tmp = (FT)1./(pc-(FT)pci)*(pc-(FT)pci);
	  res += tmp * f[index(pri,pci)]; w += tmp;
	}
	if (f[index(pri,pci+1)] != outside) // two points involved
	{
	  tmp = (FT)1./(pc-(FT)(pci+1))*(pc-(FT)(pci+1));
	  return ((res + tmp * f[index(pri,pci+1)]) / (w + tmp));
	}
	return f[index(pri,pci)]; // both are outside
      }
      if (pc-(FT)pci == 0) // same column, but not same row
      {
	if (f[index(pri,pci)] != outside)
	{
	  tmp = (FT)1./(pc-(FT)pri)*(pc-(FT)pri);
	  res += tmp * f[index(pri,pci)]; w += tmp;
	}
	if (f[index(pri+1,pci)] != outside) // two points involved
	{
	  tmp = (FT)1./(pc-(FT)(pri+1))*(pc-(FT)(pri+1));
	  return ((res + tmp * f[index(pri+1,pci)]) / (w + tmp));
	}
	return f[index(pri,pci)]; // both are outside
      }
      // --- main case -----------------------------
      IT r3 = (IT)pr, c2 = (IT)pc;
      IT r1 = r3+1, c1 = c2+1;

      IT i1 = index(r1,c1), i2 = index(r1,c2), i3 = index(r3,c2), i4 = index(r3,c1);
      if (f[i1] != outside)
      {
	tmp = (FT)1./((pr-(FT)r1)*(pr-(FT)r1)+(pc-(FT)c1)*(pc-(FT)c1));
	res += tmp * f[i1]; w += tmp;
      }
      if (f[i2]  != outside)
      {
	tmp = (FT)1./((pr-(FT)r1)*(pr-(FT)r1)+(pc-(FT)c2)*(pc-(FT)c2));
	res += tmp * f[i2]; w += tmp;
      }
      if (f[i3]  != outside)
      {
	tmp = (FT)1./((pr-(FT)r3)*(pr-(FT)r3)+(pc-(FT)c2)*(pc-(FT)c2));
	res += tmp * f[i3]; w += tmp;
      }
      if (f[i4]  != outside)
      {
	tmp = (FT)1./((pr-(FT)r3)*(pr-(FT)r3)+(pc-(FT)c1)*(pc-(FT)c1));
	res += tmp * f[i4]; w += tmp;
      }
      return (res / w);  // TODO : both outside
    }
    // ------------------------------------------------------------
    /**
     * @brief Back propagation from a given point (floatting-point coordinate), knowing the gradient of a function, by solving an ODE by Euler method.
     * @param[in] r Initial point row coordinate.
     * @param[in] c Initial point column coordinate.
     * @param[in] Dr Gradient on rows.
     * @param[in] Dc Gradient on columns.
     * @param[in,out] path Points of the solution of the ODE.
     * @param[in] dt Step size (0.6 by default).
     */
    template <typename FT, class Container>
    bool odeEuler(FT r, FT c, const FT *Dr, const FT *Dc, Container &path, bool terminal_nodes = true, FT dt = 0.6)
    {
      IT p = index((IT)(r+(FT)0.5),(IT)(c+(FT)0.5)), ri, ci, mx = size()/10;
      if (Dr[p] == 0 && Dc[p] == 0) return false;
      FT outside = std::numeric_limits<FT>::max();
      FT dr, dc, dtn;
      if (terminal_nodes) path.push_back(std::make_pair(r,c));
      for (IT brk = 0; brk < mx; brk++)
      {
        dr = nnInterpolation(r,c,Dr,outside);
        dc = nnInterpolation(r,c,Dc,outside);
        dtn = dt / std::sqrt(dr*dr + dc*dc);
        r = r - dtn * dr; c = c - dtn * dc;
        if (!(r>=0 && r<_rows && c>=0 && c<_cols))
            return false;

        path.push_back(std::make_pair(r,c));
        ri = (IT)(r+(FT)0.5);
        ci = (IT)(c+(FT)0.5);
        p = index(ri,ci);
        if (Dr[p] == 0 && Dc[p] == 0)
        {
            if (terminal_nodes) path.push_back(std::make_pair((FT)ri,(FT)ci));
            return true;
        }
      }
      return false;
    }
    // ------------------------------------------------------------
    template <typename FT>
    IT discreteMinimalPath(IT idx, const FT *Dr, const FT *Dc, const IT *V, IT mxSteps = 5000) const
    {
      unsigned short k;
      IT r, c, p = idx, pp, vold = -1;
      point(idx,r,c);
      if (Dr[p] == 0.0 && Dc[p] == 0.0) { return V[idx]; }
      FT rf = (FT)r, cf = (FT)c, dr, dc, dtn, dt = (FT)2.0, rff, cff;
      dt = std::sqrt(dt);
      for (IT brk = 0; brk < mxSteps; brk++)
      {
	dr = Dr[p]; dc = Dc[p];
	dtn = std::sqrt(dr*dr + dc*dc);
	rff = rf - dr / dtn;
	cff = cf - dc / dtn;
	r = row(rff); c = col(cff);
	pp = index(r,c);
	if (p == pp)
	{
	  rff = rf - dr * dt / dtn;
	  cff = cf - dc * dt / dtn;
	  r = row(rff); c = col(cff);
	  pp = index(r,c);
	  if (p == pp) { std::cerr << "Stop path\n"; return V[pp]; }
	}
	rf = rff; cf = cff;
	p = pp;
	if (Dr[p] == 0.0 && Dc[p] == 0.0) return V[p];
	//
	k = 0;
	for (; k < 8; k++)
	{
	  pp = neighbor8(p,k);
	  if (V[pp] == 0) continue;
	  if (V[pp] != V[p]) break;
	}
	if (k == 8) { return V[p]; }
	//
	if (vold < 0) { vold = V[p]; }
	else
	  if (V[p] != vold) { vold = V[p]; }
      }
      return vold;
    }
    // ------------------------------------------------------------
    template <typename FT, class Container1, class Container2>
    bool discreteMinimalPath(IT r, IT c, const FT *Dr, const FT *Dc, Container1 &path, Container2 &dpath,
			     bool terminal_nodes = true)
    {
      IT p = index(r,c), pp, mx = size()/10;
      if (Dr[p] == 0.0 && Dc[p] == 0.0) return false;
      FT rf = r, cf = c, dr, dc, dtn, dt = 2.0, rff, cff;
      dt = std::sqrt(dt);
      if (terminal_nodes) { path.push_back(std::make_pair(rf,cf)); dpath.push_back(p); }
      for (IT brk = 0; brk < mx; brk++)
      {
	dr = Dr[p]; dc = Dc[p];
	dtn = std::sqrt(dr*dr + dc*dc);
	rff = rf - dr / dtn;
	cff = cf - dc / dtn;
	r = row(rff); c = col(cff);
	pp = index(r,c);
	if (p == pp)
	{
	  rff = rf - dr * dt / dtn;
	  cff = cf - dc * dt / dtn;
	  r = row(rff); c = col(cff);
	  pp = index(r,c);
	  if (p == pp) { std::cerr << "Stop path\n"; return false; }
	}
	rf = rff; cf = cff;
	p = pp;
	if (Dr[p] == 0.0 && Dc[p] == 0.0)
	{
	  if (terminal_nodes)
	  {
	    path.push_back(std::make_pair(r,c));
	    dpath.push_back(p);
	  }
	  return true;
	}
	else { path.push_back(std::make_pair(rf,cf)); dpath.push_back(p); }
      }
      return false;
    }
    // ------------------------------------------------------------
    template <typename FT, class Container>
    bool discreteMinimalPath(IT r, IT c, const FT *Dr, const FT *Dc, Container &path, bool terminal_nodes = true)
    {
      IT p = index(r,c), pp, mx = size()/10;
      if (Dr[p] == 0.0 && Dc[p] == 0.0) return false;
      FT rf = r, cf = c, dr, dc, dtn, dt = 2.0, rff, cff;
      dt = std::sqrt(dt);
      if (terminal_nodes) path.push_back(p);
      for (IT brk = 0; brk < mx; brk++)
      {
	dr = Dr[p]; dc = Dc[p];
	dtn = std::sqrt(dr*dr + dc*dc);
	rff = rf - dr / dtn;
	cff = cf - dc / dtn;
	r = row(rff); c = col(cff);
	pp = index(r,c);
	if (p == pp)
	{
	  rff = rf - dr * dt / dtn;
	  cff = cf - dc * dt / dtn;
	  r = row(rff); c = col(cff);
	  pp = index(r,c);
	  if (p == pp) { std::cerr << "Stop path\n"; return false; }
	}
	rf = rff; cf = cff;
	p = pp;
	if (Dr[p] == 0.0 && Dc[p] == 0.0)
	{
	  if (terminal_nodes) path.push_back(p);
	  return true;
	}
	else path.push_back(p);
      }
      return false;
    }
    // ------------------------------------------------------------
    /**
     * @brief Reserve memory for a function on grid nodes.
     * @return A pointer to the new function on grid nodes.
     * @attention A function created with this method must be deallocated with the method free.
     * @par Example
     * @code
     double *f = g.reserveNodeFct<double>();
     // do something with f
     g.free(f);
     * @endcode
     */
    template <typename FT>
    inline FT* reserveNodeFct() { FT *f = new FT[_size]; return f; }
    // ------------------------------------------------------------
    /**
     * @brief Reserve memory for a function on grid nodes.
     * @param[in,out] f Function (array).
     * * @attention A function created with this method must be deallocated with the method free.
     * @par Example
     * @code
     double *f = NULL;
     g.reserve(f);
     // do something with f
     g.free(f);
     * @endcode
     */
    template <typename FT>
    inline void reserve(FT* &f) { f = new FT[_size]; }
    // ------------------------------------------------------------
    /**
     * @brief Free memory associated to a function on grid nodes (deallocation).
     * @param[in,out] f Function.
     */
    template <typename FT>
    inline void free(FT* &f) { if (f) { delete[] f; f = NULL; } }
    // ------------------------------------------------------------
    /**
     * @brief Copy a source function to a target function, both defined on grid nodes.
     * @param[in] f_source Source function.
     * @param[in] f_dest Target function.
     * @attention Both f_source and f_dest must be allocated arrays, for instance withmethod reserve().
     */
    template <typename FT, typename FT2>
    inline void copy(FT* f_source, FT2* f_dest) const
    {
      for (IT n = 0; n < _size; n++) f_dest[n] = (FT2)f_source[n];
    }
    // ------------------------------------------------------------
    /**
     * @brief Copy a source function to a target function, both defined on grid nodes.
     * @param[in] f_source Source function.
     * @param[in] f_dest Target function.
     * @attention Both f_source and f_dest must be allocated arrays, for instance withmethod reserve().
     */
    template <typename FT>
    inline void copy(FT* f_source, FT* f_dest) const { std::memcpy(f_dest,f_source,_size*sizeof(FT)); }
    // ------------------------------------------------------------
    /**
     * @brief Set all the elements of a function on grid nodes to a given value (for all i=0,...,size()-1, f[i]=val).
     * @param[in,out] f Function on grid nodes.
     * @param[in] val Value to set.
     */
    template <typename FT>
    inline void fill(FT* f, FT val) const
    {
      for (IT i = 0; i < _size; i++) f[i] = val;
    }
    // -------------------------------------------------------------------
    /**
     * @brief Set the value of a function (array having the size of the grid) at boundary vertices of the grid.
     * @param[in,out] f Function on the grid nodes.
     * @param[in] val Value to set at boundary grid nodes.
     */
    template <typename ET>
    void fillBoundary(ET *f, ET val) const
    {
      IT i = 0, j = _size - _cols;
      for (; i < _cols; i++) f[i] = f[i+j] = val;
      for (i = _cols; i < j; i += _cols) f[i] = f[i+_cols-1] = val;
    }
    // -------------------------------------------------------------------
    /**
     * @brief Set the value of a function (array having the size of the grid) at boundary vertices of the grid.
     * @param[in,out] f Function on the grid nodes.
     * @param[in] val Value to set at boundary grid nodes.
     */
    template <typename ET>
    void fillInterior(ET *f, ET val) const
    {
      IT pos;
      for_interior_nodes(pos,IT) { f[pos] = val; }
    }
    // -------------------------------------------------------------------
    /**
     * @brief Replace values of a function (array having the size of the grid) corresponding to a given value.
     * @param[in,out] f Function on the grid nodes.
     * @param[in] val Value to replace.
     * @param[in] new_val Value to set.
     */
    template <typename ET>
    void replace(ET *f, const ET &val, ET new_val) const
    { for (IT i = 0; i < _size; i++) if (f[i] == val) f[i] = new_val; }
    // -------------------------------------------------------------------
    /**
     * @brief Replace values of a function (array having the size of the grid) corresponding to a given value.
     * @param[in,out] f Function on the grid nodes.
     * @param[in] val Value to replace.
     * @param[in] new_val Value to set.
     */
    template <typename ET>
    inline ET difference(const IT &idx, const unsigned short &k, const unsigned short &rad, ET *f) const
    { return (f[idx]-f[neighbor(idx,rad,k)]); }
    // -------------------------------------------------------------------
    template <typename ET, class WeightClass>
    void gradientNorm(ET *f, const WeightClass &W, ET *GN, unsigned short adjacency = 8) const
    {
      unsigned short rad = nodeNeighborhoodIndex(adjacency), k;
      IT i;
      ET sum, diff;
      // interior nodes
      for (IT n = firstInteriorBlockBegin(), stop = firstInteriorBlockEnd();
	   stop < lastInteriorBlockEnd(); stop = nextInteriorBlockEnd(stop), n = nextInteriorBlockBegin(n))
	for (i = n; i < stop; i++)
	{
	  sum = 0;
	  for (k = 0; k < adjacency; k++) { diff = difference(i,k,rad,f); sum += W(i,k)*diff*diff; }
	  GN[i] = std::sqrt(sum/*/4.0*/);
	}
      // boundary nodes

    }
    // -------------------------------------------------------------------
    /**
     * @brief Norm of the gradient of a function defined on the grid nodes.
     * @param[in] f Function.
     * @param[out] GN Gradient norm of f.
     * @param[in,optional] direction Difference operator: 0=centered (default), 1=forward, 2=backward.
     */
    template <typename ET>
    void gradientNorm(ET *f, ET *GN, unsigned short direction = 0) const
    {
      IT i;
      ET sum, diff;
      // interior nodes
      switch (direction)
      {
      case 1: // forward differences
      {
	for (IT n = firstInteriorBlockBegin(), stop = firstInteriorBlockEnd();
	     stop < lastInteriorBlockEnd(); stop = nextInteriorBlockEnd(stop), n = nextInteriorBlockBegin(n))
	  for (i = n; i < stop; i++)
	  {
	    diff = f[i] - f[neighbor4(i,0)];
	    sum = diff * diff;
	    diff = f[i] - f[neighbor4(i,1)];
	    GN[i] = std::sqrt(sum+(diff*diff));
	  }
      }
      break;
      case 2: // backward differences
      {
	for (IT n = firstInteriorBlockBegin(), stop = firstInteriorBlockEnd();
	     stop < lastInteriorBlockEnd(); stop = nextInteriorBlockEnd(stop), n = nextInteriorBlockBegin(n))
	  for (i = n; i < stop; i++)
	  {
	    diff = f[i] - f[neighbor4(i,2)];
	    sum = diff * diff;
	    diff = f[i] - f[neighbor4(i,3)];
	    GN[i] = std::sqrt(sum + (diff * diff));
	  }
      }
      break;
      case 3: // forward+backward differences
      {
	for (IT n = firstInteriorBlockBegin(), stop = firstInteriorBlockEnd();
	     stop < lastInteriorBlockEnd(); stop = nextInteriorBlockEnd(stop), n = nextInteriorBlockBegin(n))
	  for (i = n; i < stop; i++)
	  {
	    diff = f[i] - f[neighbor4(i,2)];
	    sum = diff * diff;
	    diff = f[i] - f[neighbor4(i,3)];
	    sum += diff * diff;
	    diff = f[i] - f[neighbor4(i,0)];
	    sum += diff * diff;
	    diff = f[i] - f[neighbor4(i,1)];
	    GN[i] = std::sqrt((sum+(diff*diff)));
	  }
      }
      break;
      default: // centered differences
      {
	for (IT n = firstInteriorBlockBegin(), stop = firstInteriorBlockEnd();
	     stop < lastInteriorBlockEnd(); stop = nextInteriorBlockEnd(stop), n = nextInteriorBlockBegin(n))
	  for (i = n; i < stop; i++)
	  {
	    diff = f[neighbor4(i,2)] - f[neighbor4(i,0)];
	    sum = diff * diff;
	    diff = f[neighbor4(i,3)] - f[neighbor4(i,1)];
	    GN[i] = std::sqrt((sum + (diff * diff)) / 4.0);
	  }
      }
      } // end switch
    }
    // ------------------------------------------------------------
    /*template <typename FT>
    FT meanOnBall(const IndexType &i, FT *f, FT *fr, FT *fc, const IndexType &delta, const IndexType &rad)
    {
      IndexType ri, ci, rj, cj;
      point(i,ri,ci);
      IndexType sc = delta+1+rad;
      IndexType rnormal = fr[i], cnormal = fc[j];
      FT r = ri + rnormal, c = ci + cnormal, sum = 0.0, div = 0.0;
      IndexType rr, cc;
      for (FT dr = -rad, dc = 0; dr <= rad; dr+=1.0)
	for (dc = -rad; dc <= rad; dc+=1.0)
	{
	  rr = (IndexType)(r+dr); cc = (IndexType)(c+dc);
	  if (_G.onBoundary(_G.index(rr,cc))) continue;
	  if (((FT)rr-r)*((FT)rr-r)+((FT)cc-c)*((FT)cc-c) <= (FT)(rad*rad))
	  {
	    sum += fct[_G.index(rr,cc)];
	    div += 1.0;
	    //fct[_G.index(rr,cc)] = 255;
	  }
	}
      return (div > 0.0 ? (sum/div) : -1.0);
      }*/
    // ------------------------------------------------------------

  };
  // ===================================================
  // static members declaration
  // ===================================================
  template <typename IT> const unsigned short Grid<IT>::_opposite4[4] = { 2, 3, 0, 1 };
  // ------------------------------------------------------------
  template <typename IT> const unsigned short Grid<IT>::_opposite8[8] = { 4, 5, 6, 7, 0, 1, 2, 3 };
  // ------------------------------------------------------------
  template <typename IT> const unsigned short* Grid<IT>::_opposite[2] = { Grid<IT>::_opposite4, Grid<IT>::_opposite8};
  // ===================================================
  // ------------------------------------------------------------
  /**
   * @class Grid4.
   * @brief Grid structure in 4-adjacency.
   */
  template <typename IT = long>
  class Grid4 : public Grid<IT>
  {
  public:
    typedef IT IndexType;
    typedef IT Adj4;
    unsigned short adjacency;
    /**
     * @brief Constructor.
     */
    Grid4(IT r, IT c) : Grid<IT>(r,c), adjacency(4) { }
    /**
     * @brief Destructor.
     */
    ~Grid4() { }
    // ------------------------------------------------------------
    /**
     * @brief Index of the k-st neighbor of a node in 4-adjacency.
     * @param[in] idx Index of a node.
     * @param[in] k Neighbor of idx, which must be in [0,3].
     * @return Index of the k-st neighbor of a node in 4-adjacency.
     */
    inline IT neighbor(const IT &idx, const unsigned short &k) const { return Grid<IT>::neighbor4(idx,k); }
    // ------------------------------------------------------------
    inline const unsigned short& opposite(const unsigned short &k) const
    { return Grid<IT>::opposite4(k); }
  };
  // ------------------------------------------------------------
  /**
   * @class Grid8.
   * @brief Grid structure in 8-adjacency.
   */
  template <typename IT>
  class Grid8 : public Grid<IT>
  {
  public:
    typedef IT IndexType;
    typedef IT Adj8;
    unsigned short adjacency;
    /**
     * @brief Constructor.
     */
    Grid8(IT r, IT c) : Grid<IT>(r,c), adjacency(8) { }
    /**
     * @brief Destructor.
     */
    ~Grid8() { }
    // ------------------------------------------------------------
    /**
     * @brief Index of the k-st neighbor of a node in 8-adjacency.
     * @param[in] idx Index of a node.
     * @param[in] k Neighbor of idx, which must be in [0..7].
     * @return Index of the k-st neighbor of a node in 8-adjacency.
     */
    inline IT neighbor(const IT &idx, const unsigned short &k) const { return Grid<IT>::neighbor8(idx,k); }
    // ------------------------------------------------------------
    inline const unsigned short& opposite(const unsigned short &k) const
    { return Grid<IT>::opposite8(k); }

  };
  // ===================================================
  /**
   * @class GridWeight
   * @brief Data structure to store weights between nodes of the grid (on edges). This encodes a weighted directed graph.
   * @param FT Type of the weight (double by default).
   * @param IT Type of array indices (long by default).
   * @par Example
   * @code
   // initialization
   unsigned short adjacency = 8, k;         // choose adjacency (4 or 8)
   RectGrid<long> g(250,250);               // declare a grid
   unsigned short radius = g.nodeNeighborhoodIndex(adjacency);
   GridWeight<double,long> gw(g,adjacency); // declare a weight grid

   // assign weight values to interior nodes
   long i;
   for_grid_interior(g,i,long)
     for (k = 0; k < adjacency; k++)
       gw(i,k) = someDistanceOrSimilarity(i,g.neighbor(i,radius,k));

   // use the weights
   double sum = 0.0;
   for (long i = 0; i < g.size(); i++)
      for (k = 0; k < adjacency; k++)
         sum += gw(i,k);
   * @endcode
   */
  template <typename FT, class GridClass>
  class GridWeight
  {
  public:
    typedef typename GridClass::IndexType IndexType;
  protected:
    // ------------------------------------------------------------
    /**
     * @brief Array of weights such that _W[k][idx] provides the weight from node idx and its k neighbor in 4 or 8 adjacency.
     */
    FT **_W;
    // ------------------------------------------------------------
    /**
     * @brief Number of neighbors of each grid node (4 or 8 adjacency).
     */
    unsigned short _adjacency;
    // ------------------------------------------------------------
    /**
     * @brief Grid on which the weights are defined.
     */
    const GridClass &_G;
    // ------------------------------------------------------------
    /**
     * @brief Reserve memory for the weights.
     */
    void reserve()
    {
      _W = new FT*[_adjacency];
      for (unsigned short k = 0; k < _adjacency; k++)
      {
	_W[k] = new FT[_G.size()];
	_G.fill(_W[k],(FT)0.0);
      }
    }
    // ------------------------------------------------------------
  public:
    typedef FT WeightType;
    // ------------------------------------------------------------
    /**
     * @brief Constructor.
     * @param[in] g Grid on which the weights are defined.
     * @param[in] adjacency Adjacency of the grid graph (4 or 8).
     */
    GridWeight(const GridClass &g)
      : _W(NULL), _adjacency(g.adjacency), _G(g) { reserve(); }
    // ------------------------------------------------------------
    /**
     * @brief Destructor.
     */
    ~GridWeight()
    {
      if (_W)
      {
	for (unsigned short k = 0; k < _adjacency; k++) if (_W[k]) delete[] _W[k];
	delete[] _W;
      }
    }
    // ------------------------------------------------------------
    /**
     * @brief Access to the weight from the node index to its k neighbor.
     * @param[in] index Index of the node.
     * @param[in] k Index of its neighbor, which must be in [0..3] in 4-adjacency, and in [0..7] in 8-adjacency.
     * @return The value of the weight.
     */
    const FT& operator()(const IndexType &index, const unsigned short &k) const
    { return _W[k][index]; }
    // ------------------------------------------------------------
    /**
     * @brief Modify the weight value from node index to its k neighor.
     * @param[in] index Index of the node.
     * @param[in] k Index of its neighbor, which must be in [0..3] in 4-adjacency, and in [0..7] in 8-adjacency.
     * @return The value of the weight.
     */
    FT& operator()(const IndexType &index, const unsigned short &k) { return _W[k][index]; }
    // ------------------------------------------------------------
    void differences(FT *f, FT mult = 2.0)
    {
      IndexType i, j;
      unsigned short k;
      for_grid_interior(_G,i,IndexType)
      {
	for (k = 0; k < _adjacency; k++)
	{
	  j = _G.neighbor(i,k);
	  if (f[i]-f[j] > 0) _W[k][i] = f[j]/mult;
	  else
	    if (f[i]-f[j] < 0) _W[k][i] = f[j]*mult;
	    else _W[k][i] = f[j];
	}
      }
    }
    // ------------------------------------------------------------
    void symmetricDifferences(FT *f, FT regu = 1.0)
    {
      IndexType i, j;
      unsigned short k;
      for_grid_interior(_G,i,IndexType)
      {
	for (k = 0; k < _adjacency; k++)
	{
	  j = _G.neighbor(i,k);
	  _W[k][i] = _W[_G.opposite(k)][j] = std::abs(f[i]-f[j]) + regu;
	}
      }
    }
    // ------------------------------------------------------------
    void potential(FT *f)
    {
      IndexType i;
      unsigned short k;
      for_grid_interior(_G,i,IndexType)
      {
	for (k = 0; k < _adjacency; k++)
	{
	  _W[k][i] = f[_G.neighbor(i,k)];
	}
      }
    }
  };
  // ===================================================
  template <typename FT, class GridClass>
  class GridWeightExp : public GridWeight<FT,GridClass>
  {
  public:
    typedef typename GridWeight<FT,GridClass>::IndexType IndexType;
  protected:
    using GridWeight<FT,GridClass>::_G;
    using GridWeight<FT,GridClass>::_adjacency;
    using GridWeight<FT,GridClass>::_W;
    FT _sigma;
  public:
    GridWeightExp(const GridClass &g) : GridWeight<FT,GridClass>(g) { }
    ~GridWeightExp() { }
    void assignWeights(FT *f, FT sigma)
    {
      _sigma = 2.0*sigma*sigma;
      IndexType i,j;
      unsigned short k, rad = _G.nodeNeighborhoodIndex(_adjacency);
      FT diff, w;
      for_grid_interior(_G,i,IndexType)
      {
	for (k = 0; k < _adjacency; k++)
	{
	  j = _G.neighbor(i,rad,k);
	  diff = f[i] - f[j];
	  w = std::exp(-((diff*diff)/_sigma));
	  _W[_G.opposite(rad,k)][j] = w;
	  _W[k][i] = w;
	}
      }
    }
  };
  // ------------------------------------------------------------
  // ===================================================
  template <typename FT, class GridClass>
  class GridWeightBall : public GridWeight<FT,GridClass>
  {
  public:
    typedef typename GridWeight<FT,GridClass>::IndexType IndexType;
  protected:
    using GridWeight<FT,GridClass>::_G;
    using GridWeight<FT,GridClass>::_adjacency;
    using GridWeight<FT,GridClass>::_W;
    FT meanOnBall(const IndexType &i, const IndexType &j, FT *fct, const IndexType &delta, const IndexType &rad)
    {
      IndexType ri, ci, rj, cj;
      _G.point(i,ri,ci);
      _G.point(j,rj,cj);
      IndexType sc = delta+1+rad;
      IndexType rnormal = sc*(ci-cj), cnormal = sc*(rj-ri);
      FT r = ((ri+rj)/2.0) + rnormal, c = ((ci+cj)/2.0) + cnormal, sum = 0.0, div = 0.0;
      IndexType rr, cc;
      for (FT dr = -rad, dc = 0; dr <= rad; dr+=1.0)
	for (dc = -rad; dc <= rad; dc+=1.0)
	{
	  rr = (IndexType)(r+dr); cc = (IndexType)(c+dc);
	  if (_G.onBoundary(_G.index(rr,cc))) continue;
	  if (((FT)rr-r)*((FT)rr-r)+((FT)cc-c)*((FT)cc-c) <= (FT)(rad*rad))
	  {
	    sum += fct[_G.index(rr,cc)];
	    div += 1.0;
	    //fct[_G.index(rr,cc)] = 255;
	  }
	}
      return (div > 0.0 ? (sum/div) : -1.0);
    }
  public:
    GridWeightBall(const GridClass &g, FT *fct, const IndexType &delta, const IndexType &radius)
      : GridWeight<FT,GridClass>(g)
    {
      IndexType i ,j;
      unsigned short k;
      for_grid_interior(_G,i,IndexType)
      {
	for (k = 0; k < _adjacency; k++)
	{
	  j = _G.neighbor(i,k);
	  _W[k][i] = meanOnBall(i,j,fct,delta,radius);
	  //std::cerr << k << " " << _W[k][i] << "\n";
	}
      }
    }
    GridWeightBall(const GridWeightBall<FT,GridClass> &GW)
      : GridWeight<FT,GridClass>(GW._G)
    {
      IndexType i,j;
      unsigned short k;
      for_grid_interior(_G,i,IndexType)
      {
	for (k = 0; k < _adjacency; k++)
	{
	  j = _G.neighbor(i,k);
	  _W[k][i] = GW._W[k][i];
	}
      }
    }
    void leftTurnWeights(const GridWeightBall<FT,GridClass> &GW, FT *pot,
			 const FT &val_ref_int, const FT &val_ref_out)
    {
      IndexType i,j;
      unsigned short k, rad = _G.nodeNeighborhoodIndex(_adjacency), op;
      FT wi, wj, mx = std::numeric_limits<FT>::max();
      for_grid_interior(_G,i,IndexType)
      {
	for (k = 0; k < _adjacency; k++)
	{
	  j = _G.neighbor(i,k);
	  wi = GW._W[k][i]; // left
	  op = _G.opposite(k);
	  wj = GW._W[op][j]; // right
	  if (wi < 0)
	  {
	    _W[k][i] = mx;
	    //if (wj < 0) { _W[k][i] = mx; }
	    //else _W[k][i] = (!ODD(k) ? 1. : std::sqrt(2.0))+pot[j];
	    //else { _W[k][i] = (std::abs(wj-val_ref_out)+(!ODD(k) ? 1. : std::sqrt(2.0)));/*+pot[j])/2.0;*/}
	  }
	  else // wi >= 0
	  {
	    if (wj < 0) //_W[k][i] = (!ODD(k) ? 1. : std::sqrt(2.0))+pot[j];
	      _W[k][i] = mx;//(std::abs(wi-val_ref_int)+(!ODD(k) ? 1. : std::sqrt(2.0)));//+pot[j])/2.0;
	    else //_W[k][i] = (!ODD(k) ? pot[j]+0.001 : (pot[j]+0.001)*std::sqrt(2.0));
	      //_W[k][i] =  ((std::abs(wi-val_ref_int)+std::abs(wj-val_ref_out))/2.0*(!ODD(k) ? pot[j] : std::sqrt(2.0)*pot[j]));
	      _W[k][i] =  ((std::abs(wi-val_ref_int)+std::abs(wj-val_ref_out))/2.0*(!ODD(k) ? 1.0 : std::sqrt(2.0)))+0.001*(!ODD(k) ? 1.0 : std::sqrt(2.0));
	  }
	}
      }
    }
    void rightTurnWeights(const GridWeightBall<FT,GridClass> &GW, FT *pot,
			  const FT &val_ref_int, const FT &val_ref_out)
    {
      IndexType i,j;
      unsigned short k, rad = _G.nodeNeighborhoodIndex(_adjacency), op;
      FT wi, wj, mx = std::numeric_limits<FT>::max();
      for_grid_interior(_G,i,IndexType)
      {
	for (k = 0; k < _adjacency; k++)
	{
	  j = _G.neighbor(i,k);
	  wj = GW._W[k][i]; // left
	  op = _G.opposite(k);
	  wi = GW._W[op][j]; // right
	  if (wi < 0)
	  {
	    _W[k][i] = mx;
	    //if (wj < 0) { _W[k][i] = mx; }
	    //else _W[k][i] = (!ODD(k) ? 1. : std::sqrt(2.0))+pot[j];
	    //else { _W[k][i] = (std::abs(wj-val_ref_out)+(!ODD(k) ? 1. : std::sqrt(2.0)));/*+pot[j])/2.0;*/}
	  }
	  else // wi >= 0
	  {
	    if (wj < 0) //_W[k][i] = (!ODD(k) ? 1. : std::sqrt(2.0))+pot[j];
	      _W[k][i] = mx;//(std::abs(wi-val_ref_int)+(!ODD(k) ? 1. : std::sqrt(2.0)));//+pot[j])/2.0;
	    else //_W[k][i] = (!ODD(k) ? pot[j]+0.001 : (pot[j]+0.001)*std::sqrt(2.0));
	      //_W[k][i] =  ((std::abs(wi-val_ref_int)+std::abs(wj-val_ref_out))/2.0*(!ODD(k) ? pot[j] : std::sqrt(2.0)*pot[j]));
	      _W[k][i] =  ((std::abs(wi-val_ref_int)+std::abs(wj-val_ref_out))/2.0*(!ODD(k) ? 1.0 : std::sqrt(2.0)))+.001*(!ODD(k) ? 1.0 : std::sqrt(2.0));
	  }
	}
      }
    }

  };

}

#endif
