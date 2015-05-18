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
 * @file dijkstragrid.hh
 */
// =========================================================================
#ifndef __DIJKSTRA_GRID_HH__
#define __DIJKSTRA_GRID_HH__

#include <queue>
#include <list>
#include "minheap.hh"
#include "grids.hh"

namespace GeoComp
{
  // ------------------------------------------------------------
#define DIJKSTRA_UPDATE(newValue,cond)					\
  _st = _S[_idx_n];							\
  if (_st == S_FAR)							\
  {									\
    U[_idx_n] = U[idx] + newValue;					\
    cond;								\
    _hp->push(_idx_n);							\
    continue;								\
  }									\
  if (_st  < S_FAR)							\
  {									\
    _unew = U[idx] + newValue;						\
    if (_unew < U[_idx_n]) { U[_idx_n] = _unew; cond; _hp->update(_idx_n); } \
  }
  // ------------------------------------------------------------
#define GFT_UPDATE(newValue,cond)					\
  _st = _S[_idx_n];							\
  if (_st == S_FAR)							\
  {									\
    U[_idx_n] = newValue;						\
    cond;								\
    _hp->push(_idx_n);							\
    continue;								\
  }									\
  if (_st == S_TRIAL)							\
  {									\
    if (newValue < U[_idx_n]) { U[_idx_n] = newValue; cond; _hp->update(_idx_n); } \
  }
  // ------------------------------------------------------------
#define DIJKSTRA_UPDATE_BETWEEN(newValue,oppValue,cond)			\
  _st = _S[_idx_n];							\
  if (_st == S_FAR)							\
  {									\
    U[_idx_n] = U[idx] + newValue;					\
    cond;								\
    _hp->push(_idx_n);							\
    continue;								\
  }									\
  if (_st < S_FAR)							\
  {									\
    _unew = U[idx] + newValue;						\
    if (_unew < U[_idx_n]) { U[_idx_n] = _unew; cond; _hp->update(_idx_n); } \
    continue;								\
  }									\
  if (_st > S_BORDER && _st != _S[idx])					\
  {									\
    _uopp = U[_idx_n] + oppValue;					\
    if (_uopp == U[idx]) _E->push(new DijEdge<IT,FT>(_uopp,_idx_n,idx)); \
    else _E->push(new DijEdge<IT,FT>((U[_idx_n]+U[idx]+oppValue/2.0+newValue/2.0)/2.0,idx,_idx_n)); \
    res = false;							\
  }
  // ------------------------------------------------------------
   /**
   * @class DijEdge
   * @brief Edge between two grid nodes (used by DijkstraGrid).
   */
  // ------------------------------------------------------------
  template <typename IT = long, typename FT = double>
  struct DijEdge
  {
    /**
     * @brief Distance value associated to the edge (distance between 2 meeting fronts).
     */
    FT u_value;
    /**
     * @brief First node of the edge.
     */
    IT node1;
    /**
     * @brief Second node of the edge.
     */
    IT node2;
    /**
     * @brief Default constructor.
     */
    DijEdge() { }
    /**
     * @brief Constructor initializing members.
     */
    DijEdge(FT val, const IT &n1, const IT &n2) : u_value(val), node1(n1), node2(n2) { }
    /**
     * @brief Destructor.
     */
    ~DijEdge() { }
  };
  // ------------------------------------------------------------
   /**
   * @class DijEdgeCompare
   * @brief Comparison of two DijEdge (used by DijkstraGrid).
   */
  // ------------------------------------------------------------
  class DijEdgeCompare
  {
  public:
    /**
     * @brief Default constructor.
     */
    DijEdgeCompare() { }
    /**
     * @brief Destructor.
     */
    ~DijEdgeCompare() { }
    /**
     * @brief Comparison operator.
     * @param[in] e1 Pointer to a DijEdge.
     * @param[in] e2 Pointer to a DijEdge.
     * @return true if e1->u_value > e2->u_value, false otherwise.
     */
    template <typename IT, typename FT>
    bool operator()(DijEdge<FT,IT> *e1, DijEdge<FT,IT> *e2)
    { return (e1->u_value > e2->u_value); }
  };
  // ------------------------------------------------------------
  /**
   * @class DijkstraGrid
   * @brief Dijkstra-based tools: distance, minimal paths, partition
   * @param[in] FT Data type (signed floatting-point arithmetic).
   * @param[in] GrdClass Grid class (must inherit from Grid8 or Grid4).
   */
  // ------------------------------------------------------------
  template <typename FT = double, class GrdClass = Grid8<long> >
  class DijkstraGrid
  {
  public:
    /**
     * @brief Defines the GridClass used by Dijkstra's algorithm.
     */
    typedef GrdClass GridClass;
    /**
     * @brief Defines the index type.
     */
    typedef typename GridClass::IndexType IT;
  protected:
    // ------------------------------------------------------------
    /**
     * @brief Grid.
     */
    const GridClass &_G;
    /**
     * @brief State/mark values and partition labels.
     */
    IT *_S;
    /**
     * @brief Marker for nodes not yet in the heap.
     */
    static const IT S_FAR;
    /**
     * @brief Marker for border/outside nodes.
     */
    static const IT S_BORDER;
    /**
     * @brief Marker for nodes in the heap (in the front).
     */
    IT S_TRIAL;
    /**
     * @brief Marker for nodes already treated (represents the label of the node).
     */
    IT S_DEAD;
    /**
     * @brief Minimum heap.
     */
    MinHeap<FT,IT> *_hp;
    FT _diff;
    unsigned short _kvar;
    FT _unew;
    FT _uopp;
    IT _idx_n;
    IT _st;
    /**
     * @brief Square root of 2.
     */
    static const FT _SQRT2;
    /**
     * @brief Priority queue of DijEdge.
     */
    typedef typename std::priority_queue<DijEdge<IT,FT>*,std::vector<DijEdge<IT,FT>*>,DijEdgeCompare> EdgeQueue;
    /**
     * @brief DijEdge priority queue.
     */
    EdgeQueue *_E;
    // ------------------------------------------------------------
    template <class Grd>
    void updateNeighbors(const typename Grd::Adj4 &idx, FT *pot, FT *U)
    {
        FOREACH_4_NEIGHBOR(_kvar)
        {
            _idx_n = _G.neighbor4(idx,_kvar);
            DIJKSTRA_UPDATE(pot[_idx_n],_S[_idx_n]=-_S[idx])
        }
    }
    // ------------------------------------------------------------
    template <class Grd>
    void updateNeighbors(const typename Grd::Adj4 &idx, FT *pot, FT *U, IT *Pred)
    {
        FOREACH_4_NEIGHBOR(_kvar)
        {
            _idx_n = _G.neighbor4(idx,_kvar);
            DIJKSTRA_UPDATE(pot[_idx_n],_S[_idx_n]=-_S[idx];Pred[_idx_n]=idx)
        }
    }
    // ------------------------------------------------------------
    template <class Grd, class WeightClass>
    void updateNeighbors(const typename Grd::Adj4 &idx, typename WeightClass::WeightType *U, WeightClass &W)
    {
        FOREACH_4_NEIGHBOR(_kvar)
        {
            _idx_n = _G.neighbor4(idx,_kvar);
            DIJKSTRA_UPDATE(W(idx,_kvar),_S[_idx_n]=-_S[idx])
        }
    }

    // ------------------------------------------------------------
    template <class Grd, class WeightClass>
    void updateNeighbors(const typename Grd::Adj4 &idx, typename WeightClass::WeightType *U, IT *Pred, WeightClass &W)
    {
        FOREACH_4_NEIGHBOR(_kvar)
        {
            _idx_n = _G.neighbor4(idx,_kvar);
            DIJKSTRA_UPDATE(W(idx,_kvar),_S[_idx_n]=-_S[idx];Pred[_idx_n]=idx)
        }
    }
    // ------------------------------------------------------------
    template <class Grd, class WeightClass>
    void updateNeighbors(const typename Grd::Adj4 &idx, typename WeightClass::WeightType *pot,
			 typename WeightClass::WeightType *U, WeightClass &W)
    {
        FOREACH_4_NEIGHBOR(_kvar)
        {
            _idx_n = _G.neighbor4(idx,_kvar);
            DIJKSTRA_UPDATE(pot[_idx_n]*W(idx,_kvar),_S[_idx_n]=-_S[idx])
        }
    }

    // ------------------------------------------------------------
    template <class Grd, class WeightClass>
    void updateNeighbors(const typename Grd::Adj4 &idx, typename WeightClass::WeightType *pot,
			 typename WeightClass::WeightType *U, IT *Pred, WeightClass &W)
    {
      FOREACH_4_NEIGHBOR(_kvar)
      { _idx_n = _G.neighbor4(idx,_kvar); DIJKSTRA_UPDATE(pot[_idx_n]*W(idx,_kvar),_S[_idx_n]=-_S[idx];Pred[_idx_n]=idx) }
    }
    // ------------------------------------------------------------
    template <class Grd>
    void updateNeighbors(const typename Grd::Adj8 &idx, FT *pot, FT *U)
    {
      updateNeighbors<Grid4<IT> >(idx,pot,U);
      FOREACH_8DIAG_NEIGHBOR(_kvar) { _idx_n = _G.neighbor8(idx,_kvar); DIJKSTRA_UPDATE(_SQRT2*pot[_idx_n],_S[_idx_n]=-_S[idx]) }
    }
    // ------------------------------------------------------------
    template <class Grd>
    void updateNeighbors(const typename Grd::Adj8 &idx, FT *pot, FT *U, IT *Pred)
    {
      updateNeighbors<Grid4<IT> >(idx,pot,U,Pred);
      FOREACH_8DIAG_NEIGHBOR(_kvar)
      { _idx_n = _G.neighbor8(idx,_kvar); DIJKSTRA_UPDATE(_SQRT2*pot[_idx_n],_S[_idx_n]=-_S[idx];Pred[_idx_n]=idx) }
    }
    // ------------------------------------------------------------
    template <class Grd, class WeightClass>
    void updateNeighbors(const typename Grd::Adj8 &idx, typename WeightClass::WeightType *U, WeightClass &W)
    {
      FOREACH_8_NEIGHBOR(_kvar) { _idx_n = _G.neighbor8(idx,_kvar); DIJKSTRA_UPDATE(W(idx,_kvar),_S[_idx_n]=-_S[idx]) }
    }
    // ------------------------------------------------------------
    template <class Grd, class WeightClass>
    void updateNeighbors(const typename Grd::Adj8 &idx, typename WeightClass::WeightType *pot,
			 typename WeightClass::WeightType *U, WeightClass &W)
    {
      updateNeighbors<Grid4<IT> >(idx,pot,U,W);
      FOREACH_8DIAG_NEIGHBOR(_kvar)
      { _idx_n = _G.neighbor8(idx,_kvar); DIJKSTRA_UPDATE(_SQRT2*pot[_idx_n]*W(idx,_kvar),_S[_idx_n]=-_S[idx]) }
    }
    // ------------------------------------------------------------
    template <class Grd, class WeightClass>
    void updateNeighbors(const typename Grd::Adj8 &idx, typename WeightClass::WeightType *U, IT *Pred, WeightClass &W)
    {
      FOREACH_8_NEIGHBOR(_kvar)
      { _idx_n = _G.neighbor8(idx,_kvar); DIJKSTRA_UPDATE(W(idx,_kvar),_S[_idx_n]=-_S[idx];Pred[_idx_n]=idx) }
    }
    // ------------------------------------------------------------
    template <class Grd, class WeightClass>
    void updateNeighbors(const typename Grd::Adj8 &idx, typename WeightClass::WeightType *pot,
			 typename WeightClass::WeightType *U, IT *Pred, WeightClass &W)
    {
      updateNeighbors<Grid4<IT> >(idx,pot,U,Pred,W);
      FOREACH_8DIAG_NEIGHBOR(_kvar)
      { _idx_n = _G.neighbor8(idx,_kvar); DIJKSTRA_UPDATE(_SQRT2*pot[_idx_n]*W(idx,_kvar),_S[_idx_n]=-_S[idx];Pred[_idx_n]=idx) }
    }
    // ------------------------------------------------------------
    // GFT updates
    // ------------------------------------------------------------
    template <class Grd, class WeightClass>
    void updateNeighborsGFT(const typename Grd::Adj4 &idx, typename WeightClass::WeightType *U, IT *Pred, WeightClass &W)
    {
      FOREACH_4_NEIGHBOR(_kvar) { _idx_n = _G.neighbor4(idx,_kvar); GFT_UPDATE(std::max(W(idx,_kvar),U[idx]),_S[_idx_n]=-_S[idx];Pred[_idx_n]=idx) }
    }
    // ------------------------------------------------------------
    template <class Grd, class WeightClass>
    void updateNeighborsGFT(const typename Grd::Adj8 &idx, typename WeightClass::WeightType *U, IT *Pred, WeightClass &W)
    {
      updateNeighborsGFT<Grid4<IT> >(idx,U,Pred,W);
      FOREACH_8DIAG_NEIGHBOR(_kvar)
      { _idx_n = _G.neighbor8(idx,_kvar); GFT_UPDATE(std::max(W(idx,_kvar),U[idx]),_S[_idx_n]=-_S[idx];Pred[_idx_n]=idx) }
    }
    // ------------------------------------------------------------
    template <class Grd, class WeightClass>
    void updateNeighborsRiverbed(const typename Grd::Adj4 &idx, typename WeightClass::WeightType *U, IT *Pred, WeightClass &W)
    {
      FOREACH_4_NEIGHBOR(_kvar) { _idx_n = _G.neighbor4(idx,_kvar); GFT_UPDATE(W(idx,_kvar),_S[_idx_n]=S_TRIAL;Pred[_idx_n]=idx) }
    }
    // ------------------------------------------------------------
    template <class Grd, class WeightClass>
    void updateNeighborsRiverbed(const typename Grd::Adj8 &idx, typename WeightClass::WeightType *U, IT *Pred, WeightClass &W)
    {
        FOREACH_8_NEIGHBOR(_kvar)
        {
            _idx_n = _G.neighbor8(idx,_kvar);
            GFT_UPDATE(W(idx,_kvar),_S[_idx_n]=S_TRIAL;Pred[_idx_n]=idx)
        }
    }
    // ------------------------------------------------------------
    // Between updates
    // ------------------------------------------------------------
    template <class Grd, class WeightClass>
    bool updateNeighborsBetween(const typename Grd::Adj4 &idx, typename WeightClass::WeightType *U, WeightClass &W)
    {
      bool res = true;
      FOREACH_4_NEIGHBOR(_kvar)
      {
	_idx_n = _G.neighbor4(idx,_kvar);
	DIJKSTRA_UPDATE_BETWEEN(W(idx,_kvar),W(_idx_n,_G.opposite4(_kvar)),_S[_idx_n]=-_S[idx])
      }
      return res;
    }
    // ------------------------------------------------------------
    template <class Grd, class WeightClass>
    bool updateNeighborsBetween(const typename Grd::Adj4 &idx, typename WeightClass::WeightType *U, IT *Pred, WeightClass &W)
    {
      bool res = true;
      FOREACH_4_NEIGHBOR(_kvar)
      {
	_idx_n = _G.neighbor4(idx,_kvar);
	DIJKSTRA_UPDATE_BETWEEN(W(idx,_kvar),W(_idx_n,_G.opposite4(_kvar)),_S[_idx_n]=-_S[idx];Pred[_idx_n]=idx)
      }
      return res;
    }
    // ------------------------------------------------------------
    template <class Grd, class WeightClass>
    bool updateNeighborsBetween(const typename Grd::Adj8 &idx, typename WeightClass::WeightType *U, WeightClass &W)
    {
      bool res = updateNeighborsBetween<Grid4<IT> >(idx,U,W);
      FOREACH_8DIAG_NEIGHBOR(_kvar)
      {
	_idx_n = _G.neighbor8(idx,_kvar);
	DIJKSTRA_UPDATE_BETWEEN(_SQRT2*W(idx,_kvar),_SQRT2*W(_idx_n,_G.opposite8(_kvar)),_S[_idx_n]=-_S[idx])
      }
      return res;
    }
    // ------------------------------------------------------------
    template <class Grd, class WeightClass>
    bool updateNeighborsBetween(const typename Grd::Adj8 &idx, typename WeightClass::WeightType *U, IT *Pred, WeightClass &W)
    {
      bool res = updateNeighborsBetween<Grid4<IT> >(idx,U,Pred,W);
      FOREACH_8DIAG_NEIGHBOR(_kvar)
      {
	_idx_n = _G.neighbor8(idx,_kvar);
	DIJKSTRA_UPDATE_BETWEEN(_SQRT2*W(idx,_kvar),_SQRT2*W(_idx_n,_G.opposite8(_kvar)),_S[_idx_n]=-_S[idx];Pred[_idx_n]=idx)
      }
      return res;
    }
    // ------------------------------------------------------------
    void init(FT *pot, FT *U, IT *S)
    {
      IT r = 1, c, pos;
      if (_S) delete[] _S;
      if (S) _S = S;
      else _S = new IT[_G.size()];
      _G.fillBoundary(_S, S_BORDER);
      for_grid_interior(_G,pos,IT)
      {
	U[pos] = INFPLUS;
	_S[pos] = (pot[pos] < 0 ? S_BORDER : S_FAR);  // take into account the domain/mask
      }
      if (_hp) delete _hp;
      _hp = new MinHeap<FT,IT>(_G.size(),U);
    }
    // ------------------------------------------------------------
    void init(FT *pot, FT *U, IT *Pred, IT *S)
    {
      IT r = 1, c, pos;
      if (_S) delete[] _S;
      if (S) _S = S;
      else _S = new IT[_G.size()];
      _G.fillBoundary(_S, S_BORDER);
      for_grid_interior(_G,pos,IT)
      {
	U[pos] = INFPLUS;
	Pred[pos] = pos;
	_S[pos] = (pot[pos] < 0 ? S_BORDER : S_FAR);  // take into account the domain/mask
      }
      if (_hp) delete _hp;
      _hp = new MinHeap<FT,IT>(_G.size(), U);
    }
    // ------------------------------------------------------------
    void init(FT *U, IT *Pred, IT *S)
    {
      IT r = 1, c, pos;
      if (_S) delete[] _S;
      if (S) _S = S;
      else _S = new IT[_G.size()];
      _G.fillBoundary(_S, S_BORDER);
      for (; r < _G.rows()-1; r++)
	for (c = 1; c < _G.cols()-1; c++)
	{
	  pos = _G.index(r,c);
	  U[pos] = INFPLUS; Pred[pos] = pos;
	  _S[pos] = S_FAR;
	}
      if (_hp) delete _hp;
      _hp = new MinHeap<FT,IT>(_G.size(), U);
    }
    // ------------------------------------------------------------
    void init(FT *U, IT *S)
    {
      IT r = 1, c, pos;
      if (_S) delete[] _S;
      if (S) _S = S;
      else _S = new IT[_G.size()];
      _G.fillBoundary(_S, S_BORDER);
      for (; r < _G.rows()-1; r++)
	for (c = 1; c < _G.cols()-1; c++)
	{
	  pos = _G.index(r,c);
	  U[pos] = INFPLUS;
	  _S[pos] = S_FAR;
	}
      if (_hp) delete _hp;
      _hp = new MinHeap<FT,IT>(_G.size(), U);
    }
    // ------------------------------------------------------------
    inline void initEdgeQueue()
    {
      _E = new EdgeQueue;
    }
    // ------------------------------------------------------------
    inline void addNode(const IT &idx, FT *U)
    {
      U[idx] = 0;
      _S[idx] = S_TRIAL;
      _hp->push(idx);
    }
    // ------------------------------------------------------------
    void addNodes(IT *seeds_idx, IT nb_seeds, FT *U)
    {
      IT idx;
      for (IT r = 0; r < nb_seeds; r++)
      {
	idx = seeds_idx[r];
	U[idx] = 0;
	_S[idx] = S_TRIAL;
	_hp->push(idx);
      }
    }
    // ------------------------------------------------------------
    inline void addLabelNode(const IT &idx, FT *U)
    {
      U[idx] = 0;
      _S[idx] = S_TRIAL;
      _hp->push(idx);
      S_TRIAL--;
    }
    // ------------------------------------------------------------
  public:
    // ------------------------------------------------------------
    /**
     * @brief Max positive value of the type FT.
     */
    static const FT INFPLUS;
    // ------------------------------------------------------------
    /**
     * @brief Constructor.
     * @param[in] G Grid (must inherit from Grid8 or Grid4).
     */
    DijkstraGrid(const GridClass &G)
      : _G(G), _S(NULL), S_TRIAL(-2), S_DEAD(2), _hp(NULL), _kvar(0), _E(NULL) { }
    // ------------------------------------------------------------
    /**
     * @brief Destructor.
     */
    ~DijkstraGrid()
    {
      if (_hp) delete _hp;
      if (_S) delete[] _S;
      if (_E) { DijEdge<IT,FT> *e = NULL; while (!_E->empty()) { e = _E->top(); _E->pop(); delete e; } delete _E; }
    }
    // ------------------------------------------------------------
    /**
     * @brief Compute the distance to one source point idx, such that U[idx]=0 (potential-based version).
     * @param[in] idx Index of the source node of the grid graph satisfying U[idx]=0.
     * @param[in] pot Potential (scalar) field on the whole set of grid nodes. Nodes outside of the definition domain of U must have negative values (behaves like a mask).
     * @param[out] U Distance array on the set of nodes (arrival time).
     * @param[out] Pred Predecessor array on the set of nodes (not computed if set to NULL).
     */
    // CHECKED
    void distance(IT idx, FT *pot, FT *U, IT *Pred = NULL)
    {
      if (Pred)
      {
	init(pot,U,Pred,_S);
	addNode(idx,U);
	while (!_hp->empty()) { idx = _hp->pop(); _S[idx] = S_DEAD; updateNeighbors<GridClass>(idx,pot,U,Pred); }
      }
      else
      {
	init(pot,U,_S);
	addNode(idx,U);
	while (!_hp->empty()) { idx = _hp->pop(); _S[idx] = S_DEAD; updateNeighbors<GridClass>(idx,pot,U); }
      }
      delete[] _S; _S = NULL;
    }
    // ------------------------------------------------------------
    /**
     * @brief Compute the distance to one source point idx, such that U[idx]=0 (potential-based version).
     * @param[in] idx Index of the source node of the grid graph satisfying U[idx]=0.
     * @param[in] pot Potential (scalar) field on the whole set of grid nodes. Nodes outside of the definition domain of U must have negative values (behaves like a mask).
     * @param[out] U Distance array on the set of nodes (arrival time).
     * @param[out] Pred Predecessor array on the set of nodes (not computed if set to NULL).
     */
    // CHECKED
    void updateDensity(const IT &idx, IT *Pred, IT *Dens)
    {
      IT idxt = Pred[idx], idxp;
      do
      {
	idxp = idxt;
	Dens[idxp] += 1;
	idxt = Pred[idxp];
      } while (idxt != idxp);
    }

    void distance(IT idx, FT *pot, FT *U, IT *Dens, IT *Pred)
    {
      init(pot,U,Pred,_S);
      addNode(idx,U);
      _G.fill(Dens,0);
      while (!_hp->empty())
      {
	idx = _hp->pop();
	_S[idx] = S_DEAD;
	updateNeighbors<GridClass>(idx,pot,U,Pred);
	updateDensity(idx,Pred,Dens);
      }
      delete[] _S; _S = NULL;
    }
    // ------------------------------------------------------------
    /**
     * @brief Compute the distance to one source point idx, such that U[idx]=0 (potential-based version).
     * @param[in] point_set_idx Array of source nodes (given by their indicies).
     * @param[in] nb_points Number of source nodes.
     * @param[in] pot Potential (scalar) field on the whole set of grid nodes. Nodes outside of the definition domain of U must have negative values (behaves like a mask).
     * @param[out] U Distance array on the set of nodes (arrival time).
     * @param[out] Pred Predecessor array on the set of nodes (not computed if set to NULL).
     */
    // CHECKED
    void distance(IT *point_set_idx, IT nb_points, FT *pot, FT *U, IT *Pred = NULL)
    {
      IT idx;
      if (Pred)
      {
	init(pot,U,Pred,_S);
	addNodes(point_set_idx,nb_points,U);
	while (!_hp->empty()) { idx = _hp->pop(); _S[idx] = S_DEAD; updateNeighbors<GridClass>(idx,pot,U,Pred); }
      }
      else
      {
	init(pot,U,_S);
	addNodes(point_set_idx,nb_points,U);
	while (!_hp->empty()) { idx = _hp->pop(); _S[idx] = S_DEAD; updateNeighbors<GridClass>(idx,pot,U); }
      }
      delete[] _S; _S = NULL;
    }
    // ------------------------------------------------------------
    /**
     * @brief Compute the distance to one source point idx, such that U[idx]=0 (weight-based version).
     * @param[in] index Index of the source point on the grid.
     * @param[in] W Weight class (must inherit from GridWeight).
     * @param[out] U Distance array on the set of nodes (arrival time).
     * @param[out] Pred Predecessor array on the set of nodes (not computed if set to NULL).
     */
    template <class WeightClass>
    void distance(IT idx, WeightClass &W, FT *U, IT *Pred = NULL)
    {
      if (Pred)
      {
	init(U,Pred,_S); addNode(idx,U);
	while (!_hp->empty()) { idx = _hp->pop(); _S[idx] = -_S[idx]; updateNeighbors<GridClass,WeightClass>(idx,U,Pred,W); }
      }
      else
      {
	init(U,_S); addNode(idx,U);
	while (!_hp->empty()) { idx = _hp->pop(); _S[idx] = -_S[idx]; updateNeighbors<GridClass,WeightClass>(idx,U,W); }
      }
      delete[] _S; _S = NULL;
    }
    // ------------------------------------------------------------
    /**
     * @brief Compute the distance to one source point idx, such that U[idx]=0 (potential and weight-based version).
     * @param[in] index Index of the source point on the grid.
     * @param[in] W Weight class (must inherit from GridWeight).
     * @param[in] pot Potential array on the set of nodes.
     * @param[out] U Distance array on the set of nodes (arrival time).
     * @param[out] Pred Predecessor array on the set of nodes (not computed if set to NULL).
     */
    template <class WeightClass>
    void distance(IT idx, WeightClass &W, FT *pot, FT *U, IT *Pred = NULL)
    {
      if (Pred)
      {
	init(pot,U,Pred,_S); addNode(idx,U);
	while (!_hp->empty()) { idx = _hp->pop(); _S[idx] = S_DEAD; updateNeighbors<GridClass,WeightClass>(idx,pot,U,Pred,W); }
      }
      else
      {
	init(pot,U,_S); addNode(idx,U);
	while (!_hp->empty()) { idx = _hp->pop(); _S[idx] = S_DEAD; updateNeighbors<GridClass,WeightClass>(idx,pot,U,W); }
      }
      delete[] _S; _S = NULL;
    }
    // ------------------------------------------------------------
    /**
     * @brief Compute the distance to one source point idx, such that U[idx]=0 (potential-based version).
     * @param[in] idx Index of the source node of the grid graph satisfying U[idx]=0.
     * @param[in] idx_to Index of the target node of the grid graph satisfying. The propagation is stopped when idx_to is reached by thr front.
     * @param[in] pot Potential (scalar) field on the whole set of grid nodes. Nodes outside of the definition domain of U must have negative values (behaves like a mask).
     * @param[out] U Distance array on the set of nodes (arrival time).
     * @param[out] Pred Predecessor array on the set of nodes (not computed if set to NULL).
     */
    void distanceTo(IT idx, IT idx_to, FT *pot, FT *U, IT *Pred = NULL)
    {
      if (Pred)
      {
	init(pot,U,Pred,_S); addNode(idx,U);
	while (!_hp->empty())
	{
	  idx = _hp->pop(); _S[idx] = S_DEAD;
	  if (idx == idx_to) break;
	  updateNeighbors<GridClass>(idx,pot,U,Pred);
	}
      }
      else
      {
	init(pot,U,_S); addNode(idx,U);
	while (!_hp->empty())
	{
	  idx = _hp->pop(); _S[idx] = S_DEAD;
	  if (idx == idx_to) break;
	  updateNeighbors<GridClass>(idx,pot,U);
	}
      }
      while (!_hp->empty()) { idx = _hp->pop(); U[idx] = INFPLUS; }
      delete[] _S; _S = NULL;
    }
    // ------------------------------------------------------------
    /**
     * @brief Compute the distance between two source nodes (weight-based version). The propagation is stopped when minimal paths between the two nodes are found.
     * @param[in] idx1 Index of the 1st source node.
     * @param[in] idx2 Index of the 2nd source node.
     * @param[in] W Weight class (must inherit from GridWeight).
     * @param[out] U Distance array on the set of nodes (arrival time).
     * @param[out] cp List of DijEdge (pairs of nodes).
     * @param[out] Pred Predecessor array on the set of nodes (not computed if set to NULL).
     */
    template <class WeightClass>
    void distanceBetween(IT idx1, IT idx2, WeightClass &W, FT *U, std::list<DijEdge<IT,FT>*> &cp, IT *Pred = NULL)
    {
      DijEdge<IT,FT> *e = NULL;
      initEdgeQueue();
      FT dmin;
      if (Pred)
      {
	init(U,Pred,_S); addLabelNode(idx1,U); addLabelNode(idx2,U);
	while (!_hp->empty())
	{
	  idx1 = _hp->pop(); _S[idx1] = -_S[idx1];
	  if (!updateNeighborsBetween<GridClass>(idx1,U,Pred,W)) break;
	}
	while (!_hp->empty())
	{
	  idx1 = _hp->top(); e = _E->top();
	  if (e->u_value < U[idx1]) // minimal contact detected
	  {
	    _E->pop();
	    cp.push_back(e);  // add to list of contact nodes
	    dmin = e->u_value;
	    while (!_E->empty())
	    {
	      e = _E->top();
	      _E->pop();
	      if (e->u_value == dmin) cp.push_back(e);  // minimal contact is not unique
	      else { delete e; break; }
	    }
	    while (!_E->empty()) { e = _E->top(); _E->pop(); delete e; }
	    break;
	  }
	  _hp->pop();
	  _S[idx1] = -_S[idx1];
	  updateNeighborsBetween<GridClass>(idx1,U,Pred,W);
	}
      }
      else
      {
	init(U,_S); addLabelNode(idx1,U); addLabelNode(idx2,U);
	while (!_hp->empty())
	{
	  idx1 = _hp->pop(); _S[idx1] = -_S[idx1];
	  if (!updateNeighborsBetween<GridClass>(idx1,U,W)) break;
	}
	while (!_hp->empty())
	{
	  idx1 = _hp->top(); e = _E->top();
	  if (e->u_value < U[idx1]) // minimal contact detected
	  {
	    _E->pop();
	    cp.push_back(e);  // add to list of contact nodes
	    dmin = e->u_value;
	    while (!_E->empty())
	    {
	      e = _E->top();
	      _E->pop();
	      if (e->u_value == dmin) cp.push_back(e);  // minimal contact is not unique
	      else { delete e; break; }
	    }
	    while (!_E->empty()) { e = _E->top(); _E->pop(); delete e; }
	    break;
	  }
	  _hp->pop();
	  _S[idx1] = -_S[idx1];
	  updateNeighborsBetween<GridClass>(idx1,U,W);
	}
      }
      while (!_hp->empty()) { idx1 = _hp->pop(); U[idx1] = INFPLUS; }
      delete[] _S; _S = NULL;
    }
    // ------------------------------------------------------------
    /**
     * @brief Compute the distance to one source point idx, such that U[idx]=0 (weight-based version).
     * @param[in] index Index of the source point on the grid.
     * @param[in] W Weight class (must inherit from GridWeight).
     * @param[out] U Distance array on the set of nodes (arrival time).
     * @param[out] Pred Predecessor array on the set of nodes (not computed if set to NULL).
     */
    template <class WeightClass>
    void gft(IT idx, WeightClass &W, FT *U, IT *Pred, IT *V)
    {
      init(U,Pred,_S); addNode(idx,U);
      while (!_hp->empty()) { idx = _hp->pop(); _S[idx] = S_DEAD; updateNeighborsGFT<GridClass,WeightClass>(idx,U,Pred,W); }
      delete[] _S; _S = NULL;
    }
    // ------------------------------------------------------------
    /**
     * @brief Compute the distance to one source point idx, such that U[idx]=0 (weight-based version).
     * @param[in] index Index of the source point on the grid.
     * @param[in] W Weight class (must inherit from GridWeight).
     * @param[out] U Distance array on the set of nodes (arrival time).
     * @param[out] Pred Predecessor array on the set of nodes (not computed if set to NULL).
     */
    template <class WeightClass>
    void riverbed(IT idx, WeightClass &W, FT *U, IT *Pred)
    {
        init(U,Pred,_S); addNode(idx,U);
        while (!_hp->empty())
        {
            idx = _hp->pop();
            _S[idx] = S_DEAD;
            updateNeighborsRiverbed<GridClass,WeightClass>(idx,U,Pred,W);
        }
        delete[] _S; _S = NULL;
    }

    // ------------------------------------------------------------
    template <class WeightClass>
    void riverbed(IT *point_set_idx, IT nb_points, WeightClass &W, FT *U, IT *Pred, bool closed, std::vector<std::vector<int> > &vectPaths)
    {
        IT idx, idx_to; //, i;
        bool t;
        // std::vector<IT> path;
        std::queue<IT> dead;
        init(U,Pred,_S); // initialize for the 1st point

        // construct paths iteratively between pairs of successive points -------------
        for (IT i = 1; i < nb_points; i++)
        {
            idx_to = point_set_idx[i];
            idx = point_set_idx[i-1];
            addNode(idx,U);
            t = false;

            while (!_hp->empty())
            {
                idx = _hp->pop();
                _S[idx] = S_DEAD;
                if (idx == idx_to) { t = true; _S[idx] = S_BORDER; break; }
                dead.push(idx);
                updateNeighborsRiverbed<GridClass,WeightClass>(idx,U,Pred,W);
            }

            // clean S_TRIAL nodes ----------------------------------
            while (!_hp->empty())
            {
                idx = _hp->pop();
                _S[idx] = S_FAR;
                U[idx] = INFPLUS;
                Pred[idx] = idx;
            }

            // froze path only ----------------------------------
            if (t == false) { std::cerr << "riverbed: should not happen, node " << i << " not reached !\n"; }

           vectPaths.push_back(std::vector<IT>());
            // path.clear();
            minimalPathFT(idx_to,point_set_idx[i-1],Pred,vectPaths.back());

            while (!dead.empty())
            {
                idx = dead.front();
                dead.pop();
                if (_S[idx] == S_DEAD)
                {
                    _S[idx] = S_FAR;
                    Pred[idx] = idx;
                    U[idx] = INFPLUS;
                }
            }
        }

        // last path from nb_points-1 to 0 ----------------------------------
        if (closed)
        {
            idx_to = point_set_idx[0];
            idx = point_set_idx[nb_points-1];
            addNode(idx,U);
            _S[idx_to] = S_FAR;
            U[idx_to] = INFPLUS;
            while (!_hp->empty())
            {
                idx = _hp->pop();
                _S[idx] = S_DEAD;
                if (idx == idx_to) break;
                updateNeighborsRiverbed<GridClass,WeightClass>(idx,U,Pred,W);
            }
            while (!_hp->empty())
            {
                idx = _hp->pop();
                _S[idx] = S_FAR;
                U[idx] = INFPLUS;
                Pred[idx] = idx;
            }
            vectPaths.push_back(std::vector<IT>());
            // path.clear();
            minimalPathFT(idx_to,point_set_idx[nb_points-1],Pred,vectPaths.back());
        }
        delete[] _S; _S = NULL;
    }
    // ------------------------------------------------------------
    /**
     * @brief Compute the distance to one source point idx, such that U[idx]=0 (weight-based version).
     * @param[in] index Index of the source point on the grid.
     * @param[in] W Weight class (must inherit from GridWeight).
     * @param[out] U Distance array on the set of nodes (arrival time).
     * @param[out] Pred Predecessor array on the set of nodes (not computed if set to NULL).
     */
    template <class WeightClass>
    void liveWire(IT idx, WeightClass &W, FT *U, IT *Pred)
    {
      distance(idx,W,U,Pred);
    }
    // ------------------------------------------------------------
    template <class WeightClass>
    void liveWire(IT *point_set_idx, IT nb_points, WeightClass &W, FT *U, IT *Pred, bool closed, std::vector<std::vector<IT> > &vectPaths)
    {
        IT idx, idx_to; //, i;
        bool t;
        // std::vector<IT> path;
        std::queue<IT> dead;
        init(U,Pred,_S); // initialize for the 1st point

        // construct paths iteratively between pairs of successive points -------------
        for (IT i = 1; i < nb_points; i++)
        {
            idx_to = point_set_idx[i];
            idx = point_set_idx[i-1];
            addNode(idx,U);
            t = false;
            while (!_hp->empty())
            {
                idx = _hp->pop();
                _S[idx] = S_DEAD;
                if (idx == idx_to) { t = true; _S[idx] = S_BORDER; break; }
                dead.push(idx);
                updateNeighbors<GridClass,WeightClass>(idx,U,Pred,W);
            }

            // clean S_TRIAL nodes ----------------------------------
            while (!_hp->empty())
            {
                idx = _hp->pop();
                _S[idx] = S_FAR;
                U[idx] = INFPLUS;
                Pred[idx] = idx;
            }

            // froze path only ----------------------------------
            if (t == false) { std::cerr << "liveWire: should not happen, node " << i << " not reached !\n"; }

            vectPaths.push_back(std::vector<IT>());
            // path.clear();
            minimalPathFT(idx_to,point_set_idx[i-1],Pred,vectPaths.back());

            while (!dead.empty())
            {
                idx = dead.front();
                dead.pop();
                if (_S[idx] == S_DEAD)
                {
                    _S[idx] = S_FAR;
                    Pred[idx] = idx;
                    U[idx] = INFPLUS;
                }
            }
        }

        // last path from nb_points-1 to 0 ----------------------------------
        if (closed)
        {
            idx_to = point_set_idx[0];
            idx = point_set_idx[nb_points-1];
            addNode(idx,U);
            _S[idx_to] = S_FAR;
            U[idx_to] = INFPLUS;
            while (!_hp->empty())
            {
                idx = _hp->pop();
                _S[idx] = S_DEAD;
                if (idx == idx_to) break;
                updateNeighbors<GridClass,WeightClass>(idx,U,Pred,W);
            }

            while (!_hp->empty())
            {
                idx = _hp->pop();
                _S[idx] = S_FAR;
                U[idx] = INFPLUS;
                Pred[idx] = idx;
            }
            vectPaths.push_back(std::vector<IT>());
            // path.clear();
            minimalPathFT(idx_to,point_set_idx[nb_points-1],Pred,vectPaths.back());
        }
        delete[] _S; _S = NULL;
    }

    // ------------------------------------------------------------
    /**
     * @brief Path starting at idx_from and ending at the nearest node satisfying Pred[n]=n (U[n]=0).
     * @param[in] idx_from Index of the starting node.
     * @param[in] Pred Predecessor array on the set of nodes.
     * @param[out] path Minimal path (Container must have a method push_back(elt)).
     * @param[in] terminal_nodes Add terminal nodes to the path (by default) or not.
     */
    template <class Container>
    void minimalPathFT(IT idx_from, IT idx_to, IT *Pred, Container &path)
    {
      IT idx = Pred[idx_from];
      for (; idx != Pred[idx] && idx != idx_to; idx = Pred[idx]) { path.push_back(idx); _S[idx] = S_BORDER; }
      path.push_back(idx);
      _S[idx] = S_BORDER;
    }
    // ------------------------------------------------------------
    /**
     * @brief Path starting at idx_from and ending at the nearest node satisfying Pred[n]=n (U[n]=0).
     * @param[in] idx_from Index of the starting node.
     * @param[in] Pred Predecessor array on the set of nodes.
     * @param[out] path Minimal path (Container must have a method push_back(elt)).
     * @param[in] terminal_nodes Add terminal nodes to the path (by default) or not.
     */
    template <class Container>
    static void minimalPath(IT idx_from, IT *Pred, Container &path, bool terminal_nodes = true)
    {
      if (terminal_nodes) path.push_back(idx_from);
      IT idx = Pred[idx_from];
      for (; idx != Pred[idx] && idx != idx_from; idx = Pred[idx]) path.push_back(idx);
      if (terminal_nodes) path.push_back(idx);
    }
    // ------------------------------------------------------------
    /**
     * @brief Draw the path starting at idx_from and ending at the nearest node satisfying Pred[n]=n (U[n]=0).
     * @param[in] idx_from Index of the starting node.
     * @param[in] Pred Predecessor array on the set of nodes.
     * @param[in,out] C Array where the path is drawn (must have the size of the grid).
     * @param[in] value Value used to draw the path.
     * @param[in] terminal_nodes Add terminal nodes to the path (by default) or not.
     */
    template <typename CT>
    static void drawMinimalPath(IT idx_from, IT *Pred, CT *C, const CT &value, bool terminal_nodes = true)
    {
      if (terminal_nodes) C[idx_from] = value;
      IT idx = Pred[idx_from];
      for (; idx != Pred[idx]; idx = Pred[idx]) C[idx] = value;
      if (terminal_nodes) C[idx] = value;
    }
    // ------------------------------------------------------------
  }; // End class
  // ------------------------------------------------------------
  // Static members
  // ------------------------------------------------------------
  template <typename FT, typename GrdClass>
  const FT DijkstraGrid<FT,GrdClass >::_SQRT2 = std::sqrt((FT)2.0);
  // ------------------------------------------------------------
  template <typename FT, typename GrdClass>
  const typename GrdClass::IndexType DijkstraGrid<FT,GrdClass >::S_FAR = -1;
  // ------------------------------------------------------------
  template <typename FT, typename GrdClass>
  const typename GrdClass::IndexType DijkstraGrid<FT,GrdClass >::S_BORDER = 0;
  // ------------------------------------------------------------
  template <typename FT, typename GrdClass>
  const FT DijkstraGrid<FT,GrdClass >::INFPLUS = std::numeric_limits<FT>::max();
}
// ------------------------------------------------------------
#endif
