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
 * @file fm2d.hh
 */
// =========================================================================
#ifndef __FM2_HH__
#define __FM2_HH__

#include "minheap.hh"
#include "grids.hh"
#include "cmap.hh"
#include "polynomials.hh"
#include "proputils.hh"

namespace GeoComp
{
  // ------------------------------------------------------------

  // ------------------------------------------------------------
  template <typename IT, typename FT>
  class Saddle
  {
  public:
    FT u;   // value of the distance
    IT s1;  // 1st noe for saddle
    IT s2;  // 2nd node for saddle
    bool keep;
    Saddle() : u(std::numeric_limits<FT>::max()), s1(-1), s2(-1), keep(false) { }
    Saddle(const FT &uval, const IT &saddle1, const IT &saddle2) : u(uval), s1(saddle1), s2(saddle2), keep(false) { }
    ~Saddle() { }
  };

  // ------------------------------------------------------------
  // ------------------------------------------------------------
  template <typename IT, typename FT, class NdDt = Saddle<IT,FT>, class HEDt = Saddle<IT,FT> >
  class VoronoiInterface : protected IndexedHEGraph<NodeBaseData<IT,HalfEdgeData<IT,HEDt>,NdDt> >
  {
  public:
    // ------------------------------------------------------------
    typedef NodeBaseData<IT,HalfEdgeData<IT,HEDt>,NdDt> NodeClass;
    typedef IndexedHEGraph<NodeClass> HEGraphClass;
    typedef typename HEGraphClass::HEClass HEClass;
    typedef typename HEGraphClass::HEIterator HEIterator;
    typedef typename HEGraphClass::NodeIterator NodeIterator;
    // ------------------------------------------------------------
    typedef NdDt NodeDataClass;
    typedef HEDt HEDataClass;
    // ------------------------------------------------------------
    using HEGraphClass::incidentHalfEdge;
    using HEGraphClass::halfEdgesBegin;
    using HEGraphClass::halfEdgesEnd;
    using HEGraphClass::nodesBegin;
    using HEGraphClass::nodesEnd;
    using HEGraphClass::node;
    using HEGraphClass::degree;
    // ------------------------------------------------------------
  protected:
    // ------------------------------------------------------------
    using HEGraphClass::_nodes;
    using HEGraphClass::_hes;
    using HEGraphClass::isNode;
    using HEGraphClass::isNeighbor;
    // ------------------------------------------------------------
  public:
    // ------------------------------------------------------------
    VoronoiInterface() : HEGraphClass() { }
    // ------------------------------------------------------------
    ~VoronoiInterface() { }
    // ------------------------------------------------------------
    NodeClass* newNode(const IT &nidx, FT v, const IT &s1, const IT &s2)
    {
      NodeClass *res = HEGraphClass::newNode(nidx);
      res->data = new NodeDataClass(v,s1,s2);
      return res;
    }
    // ------------------------------------------------------------
    HEClass* newEdge(const IT &n1_idx, const IT &n2_idx, FT v, const IT &s1, const IT &s2)
    {
      HEClass *res = HEGraphClass::newEdge(n1_idx,n2_idx);
      res->data = new HEDataClass(v,s1,s2);
      res->opposite->data = res->data;
      return res;
    }
    // ------------------------------------------------------------
    HEClass* insert(const IT &i1, const IT &i2, const FT &v1, const FT &v2, const FT &v12, const IT &s1, const IT &s2)
    {
      // manage nodes
      NodeClass *n1 = isNode(i1), *n2 = isNode(i2);
      if (n1 == NULL) n1 = newNode(i1,v1,s1,s2);
      else { n1->data->u = v1; n1->data->s1 = s1; n1->data->s2 = s2; }//{ if (v1 < n1->data->u) { n1->data->u = v1; n1->data->s1 = s1; n1->data->s2 = s2; }}
      if (n2 == NULL) n2 = newNode(i2,v2,s1,s2);
      else { n2->data->u = v2; n2->data->s1 = s1; n2->data->s2 = s2; }//{ if (v2 < n2->data->u) { n2->data->u = v2; n2->data->s1 = s1; n2->data->s2 = s2; }}
      // manage edge
      HEClass *he = isNeighbor(i1,i2);
      if (he == NULL) he = newEdge(i1,i2,v12,s1,s2);
      else { he->data->u = v12; he->data->s1 = s1; he->data->s2 = s2; }//{ if (v12 < he->data->u) { he->data->u = v12; he->data->s1 = s1; he->data->s2 = s2; }}
      return he;
    }
    // ------------------------------------------------------------
    NodeClass* insert(const IT &i1, const FT &v1, const IT &s1, const IT &s2)
    {
      // manage nodes
      NodeClass *n1 = isNode(i1);
      if (n1 == NULL) n1 = newNode(i1,v1,s1,s2);
      else { n1->data->u = v1; n1->data->s1 = s1; n1->data->s2 = s2; }//{ if (v1 < n1->data->u) { n1->data->u = v1; n1->data->s1 = s1; n1->data->s2 = s2; }}
      return n1;
    }
    // ------------------------------------------------------------
    void degreeAnalysis(std::vector<IT> &vec1)
    {
      IT deg;
      for (NodeIterator nit = nodesBegin(), end = nodesEnd(); nit != end; ++nit)
      {
	deg = degree(nit->second);
	if (deg == 0) std::cerr << "Degree zero\n";
	if (deg == 1) vec1.push_back(nit->first);
	if (deg > 2) std::cerr << deg << " degree node\n";
      }
    }
    // ------------------------------------------------------------
    void reduce(NodeClass *n)
    {
      std::cerr << "reduce\n";
      HEClass *he = n->halfEdge;
      HEClass *op = NULL;
      FT diff = n->data->u - he->data->u;
      NodeClass *nop = NULL;
      FT dtmp;
      std::cerr << n->data->u << "\n";
      std::cerr << he->data->u << "\n";
      do
      {
	// look if nn is a min or a max ------------------
	dtmp = n->data->u - he->data->u;
	std::cerr << n->data->u << " " << he->data->u << "\n";
	if (diff > 0)
	{
	  if (dtmp < 0) // min is found
	  {
	    std::cerr << "min " << n->data->u << " " << n->data->keep << "\n";
	    diff = dtmp;
	  }
	}
	else
	{
	  if (diff < 0 && dtmp > 0) // max is found
	  {
	    std::cerr << "max " << n->data->u<< " " << n->data->keep  << "\n";
	    diff = dtmp;
	  }
	}
	// look if middle of is a min or max -------------
	op = he->opposite;
	nop = node(op);
	dtmp = op->data->u - nop->data->u;
	std::cerr << op->data->u << " " << nop->data->u << "\n";
	if (diff > 0)
	{
	  if (dtmp < 0) // min is found
	  {
	    std::cerr << "min " << op->data->u << " " << n->data->keep << "\n";
	    diff = dtmp;
	  }
	}
	else
	{
	  if (diff < 0 && dtmp > 0) // max is found
	  {
	    std::cerr << "max " << op->data->u << " " << n->data->keep << "\n";
	    diff = dtmp;
	  }
	}
	// goto next halfedge ----------------------------
	n = nop;
	he = op->next;
      } while (op != he);
    }
    // ------------------------------------------------------------
    void reduce()
    {
      IT deg;
      for (NodeIterator nit = nodesBegin(), end = nodesEnd(); nit != end; ++nit)
      {
	deg = degree(nit->second);
	if (deg == 1) {reduce(nit->second);
	  break;}
      }
    }
    // ------------------------------------------------------------
  };

  // ------------------------------------------------------------
  // ------------------------------------------------------------
  template <typename IT, typename FT, class HEDt = VoronoiInterface<IT,FT> >
  class DelNode : public NodeBase<IT,HalfEdgeData<IT,HEDt> >
  {
  public:
    typedef IT IdxType;
    typedef HalfEdgeData<IT,HEDt> HEClass;
    typedef NodeBase<IT,HEClass> ParentNodeClass;
    using ParentNodeClass::halfEdge;
    IT gridIdx;
    DelNode() : ParentNodeClass(), gridIdx(-1) { }
    DelNode(IT i) : ParentNodeClass(), gridIdx(i) { }
    ~DelNode() { }
  };

  // ------------------------------------------------------------
  // ------------------------------------------------------------
  template <typename IT, typename FT>
  class DelaunayGraph : protected HEGraph<DelNode<IT,FT> >
  {
  public:
    // ------------------------------------------------------------
    typedef DelNode<IT,FT> NodeClass;
    typedef HEGraph<NodeClass> HEGraphClass;
    typedef typename HEGraphClass::HEClass HEClass;
    typedef typename HEClass::HEDataClass HEDataClass;
    typedef typename HEGraphClass::HEIterator HEIterator;
    using HEGraphClass::isNeighbor;
    using HEGraphClass::incidentHalfEdge;
    using HEGraphClass::halfEdgesBegin;
    using HEGraphClass::halfEdgesEnd;
    // ------------------------------------------------------------
  protected:
    // ------------------------------------------------------------
    using HEGraphClass::_nodes;
    using HEGraphClass::_hes;
    // ------------------------------------------------------------
  public:
    // ------------------------------------------------------------
    DelaunayGraph() : HEGraphClass() { }
    // ------------------------------------------------------------
    ~DelaunayGraph() { }
    // ------------------------------------------------------------
    inline IT nbEdges() { return _hes.size()/2; }
    // ------------------------------------------------------------
    NodeClass* newNode(IT gridIdx)
    {
      NodeClass *res = new NodeClass(gridIdx);
      HEGraphClass::newNode(res);
      return res;
    }
    // ------------------------------------------------------------
    HEClass* newEdge(const IT &n1_idx, const IT &n2_idx)
    {
      HEClass *res = HEGraphClass::newEdge(n1_idx,n2_idx);
      res->data = new HEDataClass;
      res->opposite->data = res->data;
      return res;
    }
    // ------------------------------------------------------------
  };

  // ------------------------------------------------------------
  /**
   * @class FM2
   * @brief Eikonal solvers (fast marching, dijkstra), minimal paths.
   * @param[in] FT Data type (signed floatting-point arithmetic).
   * @param[in] GrdClass Grid class (must inherit from Grid8 or Grid4).
   */
  template <typename FT = double, class GrdClass = Grid8<long> >
  class FM2
  {
  public:
    /**
     * @brief Defines the GridClass used by algorithms.
     */
    typedef GrdClass GridClass;
    /**
     * @brief Defines the index type.
     */
    typedef typename GridClass::IndexType IT;
    // ------------------------------------------------------------
  protected:
    // ------------------------------------------------------------
    /**
     * @brief Grid.
     */
    const GridClass &_G;
    // ------------------------------------------------------------
    // typedef FT (FM2<FT>::*SolFct)(const FT&, const FT&, const FT&);
    // SolFct *_fct[2];
    // ------------------------------------------------------------
    /**
     * @brief Simplex neighbors in 4-adjacency.
     */
    IT *_sn4[2];
    // ------------------------------------------------------------
    /**
     * @brief Simplex neighbors in 8-adjacency.
     */
    IT *_sn8[2];
    IT *_wsn8[2];
    // ------------------------------------------------------------
    /**
     * @brief Simplex neighbors in both 4 and 8-adjacency.
     */
    IT *_sn[2][2];
    FT _label;
    IT DOT4[4][2];
    IT DOT8[4][3];
    IT DOT82[4][3];
    // ------------------------------------------------------------
    /**
     * @brief Square root of 2.
     */
    static const FT _SQRT2;
    // FT *_W[2];
    static const FT SIGNR4[4][2];
    static const FT SIGNC4[4][2];
    static const FT SIGNR8[8][2];
    static const FT SIGNC8[8][2];
    // ------------------------------------------------------------
    IT *_S;
    static const IT S_FAR;
    static const IT S_BORDER;
    static const IT S_MARK;
    static const IT S_DEAD_INSIDE;
    static const IT S_DEAD_OUTSIDE;
    IT S_TRIAL;
    IT S_DEAD;
    MinHeap<FT,IT> *_hp;
    FT _diff;
    FT _pot2;
    unsigned short _kvar;
    FT _unew;
    IT _idx_n;
    IT _idx_n1;
    IT _idx_n2;
    FT _un1;
    FT _un2;
    IT _st;
    IT _nbAcceptedInside;
    IT *_P1;
    IT *_P2;
    IndexedCMap<IT,Node<IT,FT> > *_medialgraph;
    // ------------------------------------------------------------
    void constructSimplex4Neighbors()
    {
      _sn4[0] = _sn[0][0] = new IT[4];
      _sn4[1] = _sn[0][1] = new IT[4];
      _sn4[0][3] = _sn4[0][1] = 1;
      _sn4[1][2] = _sn4[1][0] = -_G.cols();
      _sn4[0][2] = _sn4[0][0] = _G.cols();
      _sn4[1][3] = _sn4[1][1] = -1;
    }
    // ------------------------------------------------------------
    void constructSimplex8Neighbors()
    {
      _sn8[0] = _sn[1][0] = new IT[8];
      _sn8[1] = _sn[1][1] = new IT[8];
      _sn8[1][6] = _sn8[0][0] = _G.cols()-1;
      _sn8[1][2] = _sn8[1][0] = -_G.cols()-1;
      _sn8[0][7] = _sn8[0][1] = -1;
      _sn8[1][3] = _sn8[1][1] = -_G.cols();
      _sn8[1][4] = _sn8[0][2] = -_G.cols()+1;
      _sn8[0][5] = _sn8[0][3] = 1;
      _sn8[0][6] = _sn8[0][4] = _G.cols()+1;
      _sn8[1][7] = _sn8[1][5] = _G.cols();
    }
    // ------------------------------------------------------------
    void constructWSimplex8Neighbors()
    {
      _wsn8[0] = new IT[8];
      _wsn8[1] = new IT[8];
      _wsn8[1][6] = _wsn8[0][0] = 7;
      _wsn8[1][2] = _wsn8[1][0] = 1;
      _wsn8[0][7] = _wsn8[0][1] = 0;
      _wsn8[1][3] = _wsn8[1][1] = 2;
      _wsn8[1][4] = _wsn8[0][2] = 3;
      _wsn8[0][5] = _wsn8[0][3] = 4;
      _wsn8[0][6] = _wsn8[0][4] = 5;
      _wsn8[1][7] = _wsn8[1][5] = 6;
    }
    // ------------------------------------------------------------
    void constructSimplexNeighbors()
    { constructSimplex4Neighbors(); constructSimplex8Neighbors();
      // _W[0] = new FT[4]; _W[1] = new FT[8];
      // _W[0][0] = _W[0][1] = _W[0][2] = _W[0][3] = 1.0;
      // _W[1][0] = _W[1][2] = _W[1][4] = _W[1][6] = 1.0;
      // _W[1][1] = _W[1][3] = _W[1][5] = _W[1][7] = std::sqrt(2.0);
    }
    // ------------------------------------------------------------
    void constructDots()
    {
      DOT4[0][1] = DOT4[1][0] = 0;
      DOT4[0][0] = DOT4[3][1] = -_G.cols();
      DOT4[1][1] = DOT4[2][0] = -1;
      DOT4[2][1] = DOT4[3][0] = -_G.cols()-1;
      DOT8[0][0] = 0;
      DOT8[0][1] = 1;
      DOT8[0][2] = _G.cols();
      DOT8[1][0] = -1;
      DOT8[1][1] = _G.cols()-1;
      DOT8[1][2] = -2;
      DOT8[2][0] = -_G.cols()-1;
      DOT8[2][1] = -_G.cols()-2;
      DOT8[2][2] = -2*_G.cols()-1;
      DOT8[3][0] = -_G.cols();
      DOT8[3][1] = -2*_G.cols();
      DOT8[3][2] = -_G.cols()+1;
      DOT82[0][0] = 0;
      DOT82[0][1] = -_G.cols();
      DOT82[0][2] = -1;
      DOT82[1][0] = -1;
      DOT82[1][1] = 0;
      DOT82[1][2] = -_G.cols()-1;
      DOT82[2][0] = -_G.cols()-1;
      DOT82[2][1] = -1;
      DOT82[2][2] = -_G.cols();
      DOT82[3][0] = -_G.cols();
      DOT82[3][1] = -_G.cols()-1;
      DOT82[3][2] = 0;
    }

    // ------------------------------------------------------------
    /*void constructSolutions()
    {
      _fct[0] = new SolFct[4];
      _fct[0][0] = _fct[0][1] = _fct[0][2] = _fct[0][3] = &FM2<FT>::fmSolution4;
      _fct[1] = new SolFct[8];
      _fct[1][0] = _fct[1][2] = _fct[1][4] = _fct[1][6] = &FM2<FT>::fmSolution8;
      _fct[1][1] = _fct[1][3] = _fct[1][5] = _fct[1][7] = &FM2<FT>::fmSolution8diag;
      _W[0] = new FT[4]; _W[1] = new FT[8];
      _W[0][0] = _W[0][1] = _W[0][2] = _W[0][3] = 1.0;
      _W[1][0] = _W[1][2] = _W[1][4] = _W[1][6] = 1.0;
      _W[1][1] = _W[1][3] = _W[1][5] = _W[1][7] = std::sqrt(2.0);
      }*/
    // ------------------------------------------------------------
    // u1 must be a 8-neighbor, and thus u2 a 4-neighbor
    inline FT fmSolution8(const FT &pot, const FT &u1, const FT &u2)
    {
      _diff = u2 - u1;
      _diff *= _diff;
      _pot2 = pot * pot;
      return (_pot2 > 2.0 * _diff ? u2 + std::sqrt(_pot2 - _diff) : u1 + _SQRT2 * pot);
    }
    // ------------------------------------------------------------
    inline FT fmSolution8diag(const FT &pot, const FT &u1, const FT &u2)
    { return fmSolution8(pot,u2,u1); } // swap the role of u1 and u2.
    // ------------------------------------------------------------
    inline FT fmSolution4(const FT &pot, const FT &u1, const FT &u2)
    {
      _diff = u1-u2;
      return (((u1+u2)+std::sqrt(2.0*(pot*pot)-_diff*_diff))/2.0);
    }
    // ------------------------------------------------------------
    /*inline FT fmSolution(const FT &pot, const FT &u1, const FT &u2,
			 const unsigned short &radius, const unsigned short &k)
			 { return (this->*_fct[radius][k])(pot,u1,u2); }*/
    // ------------------------------------------------------------
    /*void fmUpdateNeighbors(const IT &idx, FT *pot, FT *U,
			   const unsigned short &adjacency, const unsigned short &radius)
    {
      IT st;
      IT idx_n, idx_n1, idx_n2;
      FT unew, un1, un2;
      // -- update neighbors of idx -----------------
      for (unsigned short k = 0; k < adjacency; k++)
      {
	idx_n = neighbor(idx,radius,k);
	st = _S[idx_n];
	// -- not yet treated ------------------
	if (st == S_FAR)
	{
	  _S[idx_n] = S_TRIAL;
	  U[idx_n] = U[idx] + _W[radius][k] * pot[idx_n];
	  _hp->push(idx_n);
	  continue;
	}
	// --  already treated (in the heap), try to reduce --
	if (st == S_TRIAL)
	{
	  idx_n1 = simplexNeighbor(idx_n,radius,k,0);
	  idx_n2 = simplexNeighbor(idx_n,radius,k,1);
	  unew = INFPLUS;
	  if (_S[idx_n1] == S_DEAD)
	  {
	    if (_S[idx_n2] == S_DEAD)
	    {
	      un1 = U[idx_n1];
	      un2 = U[idx_n2];
	      if (un1 < un2) unew = fmSolution(pot[idx_n],un1,U[idx],radius,k);
	      else unew = fmSolution(pot[idx_n],un2,U[idx],radius,k);
	    }
	    else unew = fmSolution(pot[idx_n],U[idx_n1],U[idx],radius,k);
	  }
	  else // idx_n1 is not S_DEAD
	  {
	    if (_S[idx_n2] == S_DEAD)
	      unew = fmSolution(pot[idx_n],U[idx_n2],U[idx],radius,k);
	    else unew = U[idx] + _W[radius][k] * pot[idx_n];
	  }
	  if (unew < U[idx_n])
	  {
	    U[idx_n] = unew;
	    _hp->update(idx_n);
	  }
	}
      }
      }*/
    // ------------------------------------------------------------
    // u1 <= u2
    bool fmSolution4(const IT &pos, const FT &pot, const FT &u1, const FT &u2,
		     FT &u, FT &dur, FT &duc, const unsigned short &k, unsigned short l)
    {
      /*if (u1 == u2)
      {
	FT disc = std::sqrt(2.0*(pot*pot));
	_unew = u1 + (disc/2.0);
	if (_unew < u)  // smaller -> update
	{
	  FT du = (_unew - u2) / (2.0*std::sqrt(0.5));
	  if (ODD(k)) { dur = du * SIGNR4[k][l]; duc = du * SIGNC4[k][l]; }
	  else { dur = du * SIGNR4[k][l]; duc = du * SIGNC4[k][l]; }
	  u = _unew;
	  _hp->update(pos);
	  return true;
	}
	return false;
	}*/
      _diff = u2 - u1;
      FT disc = std::sqrt(2.0*(pot*pot)-(_diff*_diff));
      _unew = ((u1 + u2) + disc)/2.0;
      if (_unew < u)  // smaller -> update
      {
	FT tn = 0.5 + _diff / (2.0*disc);
	FT tpos = 1.0 - tn;
	FT du = (_unew - u2 + tn*_diff)/*(tn*u1+tpos*u2))*/ / std::sqrt(tn*tn+tpos*tpos);
	if (ODD(k)) { dur = du * tpos * SIGNR4[k][l]; duc = du * tn * SIGNC4[k][l]; }
	else { dur = du * tn * SIGNR4[k][l]; duc = du * tpos * SIGNC4[k][l]; }
	u = _unew;
	_hp->update(pos);
	return true;
      }
      return false;
    }
    // ------------------------------------------------------------
    // u1 <= u2
    bool fmSolution4(const IT &pos, const FT &w1, const FT &w2, const FT &u1, const FT &u2,
		     FT &u, FT &dur, FT &duc, const unsigned short &k, unsigned short l)
    {
      _diff = u2 - u1;
      FT w12 = w1 - w2;
      FT a = 2.0*(w12*w12), b = 4.0*w2*w12, c = 2*w1*w2-(_diff*_diff), d = (_diff*_diff)-2.0*w2*w2, e = (w2*w2-(_diff*_diff))/2.0;
      FT t = quartic_newton(a,b,c,d,e,0.5,1.0e-15,20);
      FT t2 = 1.0-t;
      _unew = t*u1 + t2*u2 + (w1*t + w2*t2) * std::sqrt(t*t+t2*t2);
      // if (_unew < u)  // smaller -> update
      // {
      // 	FT du = (_unew - u2 + tn*_diff)/*(tn*u1+tpos*u2))*/ / std::sqrt(tn*tn+tpos*tpos);
      // 	if (ODD(k)) { dur = du * tpos * SIGNR4[k][l]; duc = du * tn * SIGNC4[k][l]; }
      // 	else { dur = du * tn * SIGNR4[k][l]; duc = du * tpos * SIGNC4[k][l]; }
      // 	u = _unew;
      // 	_hp->update(pos);
      // 	return true;
      // }
      return false;
    }
    // ------------------------------------------------------------
    // u1 <= u2
    bool fmSolution4(const IT &pos, const FT &pot, const FT &u1, const FT &u2,
    		     FT &u, FT &dur, FT &duc, const unsigned short &k, unsigned short l, const IT &idx1, const IT &idx2)
    {
      _diff = u2 - u1;
      FT disc = std::sqrt(2.0*(pot*pot)-(_diff*_diff));
      _unew = ((u1 + u2) + disc)/2.0;
      if (_unew < u)  // smaller -> update
      {
    	FT tn = 0.5 + _diff / (2.0*disc);
    	FT tpos = 1.0 - tn;
    	FT du = (_unew - u2 + tn*_diff) / std::sqrt(tn*tn+tpos*tpos);
    	if (ODD(k)) { dur = du * tpos * SIGNR4[k][l]; duc = du * tn * SIGNC4[k][l]; }
    	else { dur = du * tn * SIGNR4[k][l]; duc = du * tpos * SIGNC4[k][l]; }
    	u = _unew;
    	_hp->update(pos);
	if (tn == 0.0) { _P1[pos] = idx2; _P2[pos] = -1; }
	else
	  if (tn == 1.0) { _P1[pos] = idx1; _P2[pos] = -1; }
	  else { _P1[pos] = idx1; _P2[pos] = idx2; }
    	return true;
      }
      return false;
    }
    // ------------------------------------------------------------
    // u1 <= u2
    bool fmSolution4(const IT &pos, const FT &pot, const FT &u1, const FT &u2,
		     FT &u, FT &dur, FT &duc, const unsigned short &k, unsigned short l, IT &v, const IT &v1, const IT &v2)
    {
      _diff = u2 - u1;
      FT disc = std::sqrt(2.0*(pot*pot)-(_diff*_diff));
      _unew = ((u1 + u2) + disc)/2.0;
      if (_unew < u)  // smaller -> update
      {
	FT tn = 0.5 + _diff / (2.0*disc);
	FT tpos = 1.0 - tn;
	FT du = (_unew - u2 + tn*_diff)/*(tn*u1+tpos*u2))*/ / std::sqrt(tn*tn+tpos*tpos);
	if (ODD(k)) { dur = du * tpos * SIGNR4[k][l]; duc = du * tn * SIGNC4[k][l]; }
	else { dur = du * tn * SIGNR4[k][l]; duc = du * tpos * SIGNC4[k][l]; }
	u = _unew;
	_hp->update(pos);
	v = v1;
	return true;
      }
      return false;
    }
    // ------------------------------------------------------------
    // u1 <= u2
    bool fmSolution4(const IT &pos, const FT &pot, const FT &u1, const FT &u2,
		     FT &u, FT &dur, FT &duc, const unsigned short &k, unsigned short l, const IT &v1, const IT &v2, FT *DUr, FT *DUc)
    {
      _diff = u2 - u1;
      FT disc = std::sqrt(2.0*(pot*pot)-(_diff*_diff));
      _unew = ((u1 + u2) + disc)/2.0;
      if (_unew < u)  // smaller -> update
      {
	FT tn = 0.5 + _diff / (2.0*disc);
	FT tpos = 1.0 - tn;
	FT du = (_unew - u2 + tn*_diff)/*(tn*u1+tpos*u2))*/ / std::sqrt(tn*tn+tpos*tpos);
	if (ODD(k)) { dur = du * tpos * SIGNR4[k][l]; duc = du * tn * SIGNC4[k][l]; }
	else { dur = du * tn * SIGNR4[k][l]; duc = du * tpos * SIGNC4[k][l]; }
	u = _unew;
	_hp->update(pos);
	_S[pos] = -v1;
	return true;
      }
      return false;
    }
    // ------------------------------------------------------------
    // u1 <= u2, u1 8-neighbor (diagonal) to pos, u2 4-neighbor to pos
    bool fmSolution8(const IT &pos, const FT &pot, const FT &u1, const FT &u2,
		     FT &u, FT &dur, FT &duc, const unsigned short &k, unsigned short l)
    {
      _diff = u2 - u1;
      FT du;
      if (pot > _SQRT2 * _diff)
      {
	_pot2 = std::sqrt(pot*pot - _diff*_diff);
	_unew = u2 + _pot2;
	if (_unew < u)
	{
	  FT tn = _diff / _pot2;
	  FT tpos = (FT)1. - tn;
	  du = (_unew - (tn*u1+tpos*u2)) / std::sqrt((FT)1.+tn*tn);
	  if (ODD(k)) { dur = du * SIGNR8[_kvar][l]; duc = du * tn * SIGNC8[_kvar][l]; }
	  else { dur = du * tn * SIGNR8[_kvar][l]; duc = du * SIGNC8[_kvar][l]; }
	  u = _unew;
	  _hp->update(pos);
	  return true;
	}
	return false;
      }
      // else
      _unew = u1 + _SQRT2 * pot;
      if (_unew < u)
      {
	dur = pot * SIGNR8[_kvar][l]; duc = pot * SIGNC8[_kvar][l];
	u = _unew;
	_hp->update(pos);
	return true;
      }
      return false;
    }
    // ------------------------------------------------------------
    // u1 <= u2, u1 8-neighbor (diagonal) to pos, u2 4-neighbor to pos
    bool fmSolution8(const IT &pos, const FT &wj, const FT &wk, const FT &u1, const FT &u2,
		     FT &u, FT &dur, FT &duc, const unsigned short &k, unsigned short l)
    {
      _diff = u2 - u1;
      FT wjk = wj - wk, du;
      FT a = 4.0*(wjk*wjk), b = 4.0*wjk*wk, c = wk*wk + 4.0*(wjk*wjk) - _diff*_diff, d = 2.0*wjk*wk, e = wjk*wjk - _diff*_diff;
      FT t = quartic_newton(a,b,c,d,e,0.5,1.0e-15,20);
      if (t <= 0.0)
      {
	_unew = u2 + wk;
	du = wk;
	t = 0.0;
      }
      else
      {
	if (t >= 1.0)
	{
	  _unew = u1 + _SQRT2*wj;
	  du = wj;
	  t = 1.0;
	}
	else
	{
	  _unew = u2 -_diff*t + (wjk*t + wk) * std::sqrt(1.0 + t*t);
	  du = (_unew + _diff*t - u2) / std::sqrt(1.0 + t*t);
	}
      }
      if (_unew < u)
      {
	if (ODD(k)) { dur = du * SIGNR8[_kvar][l]; duc = du * t * SIGNC8[_kvar][l]; }
	else { dur = du * t * SIGNR8[_kvar][l]; duc = du * SIGNC8[_kvar][l]; }
	u = _unew;
	_hp->update(pos);
	return true;
      }
      return false;
    }
    // ------------------------------------------------------------
    // u1 <= u2, u1 8-neighbor (diagonal) to pos, u2 4-neighbor to pos
    bool fmSolution8(const IT &pos, const FT &pot, const FT &u1, const FT &u2,
		     FT &u, FT &dur, FT &duc, const unsigned short &k, unsigned short l, IT &v, const IT &v1, const IT &v2)
    {
      _diff = u2 - u1;
      FT du;
      if (pot > _SQRT2 * _diff)
      {
	_pot2 = std::sqrt(pot*pot - _diff*_diff);
	_unew = u2 + _pot2;
	if (_unew < u)
	{
	  FT tn = _diff / _pot2;
	  FT tpos = 1. - tn;
	  du = (_unew - (tn*u1+tpos*u2)) / std::sqrt(1.+tn*tn);
	  if (ODD(k)) { dur = du * SIGNR8[_kvar][l]; duc = du * tn * SIGNC8[_kvar][l]; }
	  else { dur = du * tn * SIGNR8[_kvar][l]; duc = du * SIGNC8[_kvar][l]; }
	  u = _unew;
	  _hp->update(pos);
	  v = (tn > 0.5 ? v1 : v2);
	  return true;
	}
	return false;
      }
      // else
      _unew = u1 + _SQRT2 * pot;
      if (_unew < u)
      {
	dur = pot * SIGNR8[_kvar][l]; duc = pot * SIGNC8[_kvar][l];
	u = _unew;
	_hp->update(pos);
	v = v1;
	return true;
      }
      return false;
    }
    // ------------------------------------------------------------
    // u1 <= u2, u1 8-neighbor (diagonal) to pos, u2 4-neighbor to pos
    bool fmSolution8(const IT &pos, const FT &pot, const FT &u1, const FT &u2,
		     FT &u, FT &dur, FT &duc, const unsigned short &k, unsigned short l, const IT &idx1, const IT &idx2)
    {
      _diff = u2 - u1;
      FT du;
      if (pot > _SQRT2 * _diff)
      {
	_pot2 = std::sqrt(pot*pot - _diff*_diff);
	_unew = u2 + _pot2;
	if (_unew < u)
	{
	  FT tn = _diff / _pot2;
	  FT tpos = (FT)1. - tn;
	  du = (_unew - (tn*u1+tpos*u2)) / std::sqrt((FT)1.+tn*tn);
	  if (ODD(k)) { dur = du * SIGNR8[_kvar][l]; duc = du * tn * SIGNC8[_kvar][l]; }
	  else { dur = du * tn * SIGNR8[_kvar][l]; duc = du * SIGNC8[_kvar][l]; }
	  u = _unew;
	  _hp->update(pos);
	  if (tn == 0.0) { _P1[pos] = idx2; _P2[pos] = -1; }
	  else
	    if (tn == 1.0) { _P1[pos] = idx1; _P2[pos] = -1; }
	    else { _P1[pos] = idx1; _P2[pos] = idx2; }
	  return true;
	}
	return false;
      }
      // else
      _unew = u1 + _SQRT2 * pot;
      if (_unew < u)
      {
	dur = pot * SIGNR8[_kvar][l]; duc = pot * SIGNC8[_kvar][l];
	u = _unew;
	_hp->update(pos);
	_P1[pos] = idx1; _P2[pos] = -1;
	return true;
      }
      return false;
    }
    // ------------------------------------------------------------
    // u1 <= u2, u1 8-neighbor (diagonal) to pos, u2 4-neighbor to pos
    bool fmSolution8(const IT &pos, const FT &pot, const FT &u1, const FT &u2,
		     FT &u, FT &dur, FT &duc, const unsigned short &k, unsigned short l, const IT &v1, const IT &v2, FT *DUr, FT *DUc)
    {
      _diff = u2 - u1;
      FT du;
      if (pot > _SQRT2 * _diff)
      {
	_pot2 = std::sqrt(pot*pot - _diff*_diff);
	_unew = u2 + _pot2;
	if (_unew < u)
	{
	  FT tn = _diff / _pot2;
	  FT tpos = (FT)1. - tn;
	  du = (_unew - (tn*u1+tpos*u2)) / std::sqrt((FT)1.+tn*tn);
	  if (ODD(k)) { dur = du * SIGNR8[_kvar][l]; duc = du * tn * SIGNC8[_kvar][l]; }
	  else { dur = du * tn * SIGNR8[_kvar][l]; duc = du * SIGNC8[_kvar][l]; }
	  u = _unew;
	  _hp->update(pos);
	  if (v1 == v2) _S[pos] = -v1;
	  else //_S[pos] = (tn > 0.5 ? -v1 : -v2);
	  {
	    IT ss = _G.discreteMinimalPath(pos,DUr,DUc,_S);
	    if (ss <= 0) _S[pos] = (tn > 0.5 ? -v1 : -v2);
	    else {_S[pos] = -ss; }
	  }
	  return true;
	}
	return false;
      }
      // else
      _unew = u1 + _SQRT2 * pot;
      if (_unew < u)
      {
	dur = pot * SIGNR8[_kvar][l]; duc = pot * SIGNC8[_kvar][l];
	u = _unew;
	_hp->update(pos);
	_S[pos] = -v1;
	return true;
      }
      return false;
    }
    // ------------------------------------------------------------
    /*void fmUpdateNeighbors4(const IT &idx, FT *pot, FT *U)
    {
      IT st;
      IT idx_n, idx_n1, idx_n2;
      FT unew, un1, un2;
      // -- update neighbors of idx -----------------
      for (unsigned short k = 0; k < 4; k++)
      {
	idx_n = neighbor4(idx,k);
	st = _S[idx_n];
	// -- not yet treated ------------------
	if (st == S_FAR)
	{
	  _S[idx_n] = S_TRIAL;
	  U[idx_n] = U[idx] + pot[idx_n];
	  _hp->push(idx_n);
	  continue;
	}
	// --  already treated (in the heap), try to reduce --
	if (st == S_TRIAL)
	{
	  idx_n1 = simplexNeighbor4(idx_n,k,0);
	  idx_n2 = simplexNeighbor4(idx_n,k,1);
	  unew = INFPLUS;
	  if (_S[idx_n1] == S_DEAD)
	  {
	    if (_S[idx_n2] == S_DEAD)
	    {
	      un1 = U[idx_n1]; un2 = U[idx_n2];
	      if (un1 < un2) unew = fmSolution4(pot[idx_n],un1,U[idx]);
	      else unew = fmSolution4(pot[idx_n],un2,U[idx]);
	    }
	    else unew = fmSolution4(pot[idx_n],U[idx_n1],U[idx]);
	  }
	  else // idx_n1 is not S_DEAD
	  {
	    if (_S[idx_n2] == S_DEAD) unew = fmSolution4(pot[idx_n],U[idx_n2],U[idx]);
	    else unew = U[idx] + pot[idx_n];
	  }
	  if (unew < U[idx_n]) { U[idx_n] = unew; _hp->update(idx_n); }
	}
      }
      }*/
    // ============================================================
    // update macro
    // ============================================================
#define FM_UPDATE_TRIAL(simplexFct,solFct)				\
    _idx_n1 = simplexFct(_idx_n,_kvar,0);				\
    _idx_n2 = simplexFct(_idx_n,_kvar,1);				\
    _unew = INFPLUS;							\
    if (_S[_idx_n1] > S_BORDER)						\
    {									\
      if (_S[_idx_n2] > S_BORDER)					\
      {								        \
	_un1 = U[_idx_n1]; _un2 = U[_idx_n2];				\
	if (_un1 < _un2) _unew = solFct(pot[_idx_n],_un1,U[idx]);	\
	else _unew = solFct(pot[_idx_n],_un2,U[idx]);			\
      }									\
      else _unew = solFct(pot[_idx_n],U[_idx_n1],U[idx]);		\
    }									\
    else								\
    {									\
      if (_S[_idx_n2] > S_BORDER) _unew = solFct(pot[_idx_n],U[_idx_n2],U[idx]); \
      else _unew = U[idx] + pot[_idx_n];				\
    }									\
    if (_unew < U[_idx_n]) { U[_idx_n] = _unew; _hp->update(_idx_n); }
    // ------------------------------------------------------------
#define FM_UPDATE_TRIAL_DU(simplexFct,solFct,signR,signC,k)		\
    _idx_n1 = simplexFct(_idx_n,_kvar,0);				\
    _idx_n2 = simplexFct(_idx_n,_kvar,1);				\
    _unew = INFPLUS;							\
    if (_S[_idx_n1] > S_BORDER)						\
    {									\
      if (_S[_idx_n2] > S_BORDER)					\
      {								        \
	_un1 = U[_idx_n1]; _un2 = U[_idx_n2];				\
	if (_un1 < _un2) solFct(_idx_n,pot[_idx_n],_un1,U[idx],		\
				U[_idx_n],DUr[_idx_n],DUc[_idx_n],k,0); \
	else solFct(_idx_n,pot[_idx_n],_un2,U[idx],			\
		    U[_idx_n],DUr[_idx_n],DUc[_idx_n],k,1);		\
      }									\
      else solFct(_idx_n,pot[_idx_n],U[_idx_n1],U[idx],			\
		  U[_idx_n],DUr[_idx_n],DUc[_idx_n],k,0);		\
    }									\
    else								\
    {									\
      if (_S[_idx_n2] > S_BORDER)					\
	solFct(_idx_n,pot[_idx_n],U[_idx_n2],U[idx],			\
	       U[_idx_n],DUr[_idx_n],DUc[_idx_n],k,1);			\
      else								\
      {								        \
	_unew = U[idx] + pot[_idx_n];					\
	if (_unew < U[_idx_n])						\
	{								\
	  U[_idx_n] = _unew;						\
	  if (ODD(k)) { DUr[_idx_n] = pot[_idx_n] * signR[_kvar][1]; DUc[_idx_n] = 0; } \
	  else { DUr[_idx_n] = 0; DUc[_idx_n] = pot[_idx_n] * signC[_kvar][1]; } \
	  _hp->update(_idx_n);						\
	}								\
      }									\
    }
    // ------------------------------------------------------------
#define FM_UPDATE_TRIAL_DU_V(simplexFct,solFct,signR,signC,k,V)		\
    _idx_n1 = simplexFct(_idx_n,_kvar,0);				\
    _idx_n2 = simplexFct(_idx_n,_kvar,1);				\
    _unew = INFPLUS;							\
    if (_S[_idx_n1] > S_BORDER)						\
    {									\
      if (_S[_idx_n2] > S_BORDER)					\
      {								        \
	_un1 = U[_idx_n1]; _un2 = U[_idx_n2];				\
	if (_un1 < _un2)						\
	{								\
	  solFct(_idx_n,pot[_idx_n],_un1,U[idx],			\
		 U[_idx_n],DUr[_idx_n],DUc[_idx_n],k,0,V[_idx_n],V[_idx_n1],V[idx]); \
	}								\
	else								\
	{								\
	  solFct(_idx_n,pot[_idx_n],_un2,U[idx],			\
		 U[_idx_n],DUr[_idx_n],DUc[_idx_n],k,1,V[_idx_n],V[_idx_n2],V[idx]);		\
	}								\
      }									\
      else								\
      {								\
	  solFct(_idx_n,pot[_idx_n],U[_idx_n1],U[idx],			\
		 U[_idx_n],DUr[_idx_n],DUc[_idx_n],k,0,V[_idx_n],V[_idx_n1],V[idx]);		\
      }									\
    }									\
    else								\
    {									\
      if (_S[_idx_n2] > S_BORDER)					\
      {								\
	solFct(_idx_n,pot[_idx_n],U[_idx_n2],U[idx],			\
	       U[_idx_n],DUr[_idx_n],DUc[_idx_n],k,1,V[_idx_n],V[_idx_n2],V[idx]);			\
      }									\
      else								\
      {								        \
	_unew = U[idx] + pot[_idx_n];					\
	if (_unew < U[_idx_n])						\
	{								\
	  U[_idx_n] = _unew;						\
	  V[_idx_n] = V[idx];						\
	  if (ODD(k)) { DUr[_idx_n] = pot[_idx_n] * signR[_kvar][1]; DUc[_idx_n] = 0; } \
	  else { DUr[_idx_n] = 0; DUc[_idx_n] = pot[_idx_n] * signC[_kvar][1]; } \
	  _hp->update(_idx_n);						\
	}								\
      }									\
    }
     // ------------------------------------------------------------
#define FM_UPDATE_TRIAL_DU_SV(simplexFct,solFct,signR,signC,k)		\
    _idx_n1 = simplexFct(_idx_n,_kvar,0);				\
    _idx_n2 = simplexFct(_idx_n,_kvar,1);				\
    _unew = INFPLUS;							\
    if (_S[_idx_n1] > S_BORDER)						\
    {									\
      if (_S[_idx_n2] > S_BORDER)					\
      {								        \
	_un1 = U[_idx_n1]; _un2 = U[_idx_n2];				\
	if (_un1 < _un2)						\
	{								\
	  solFct(_idx_n,pot[_idx_n],_un1,U[idx],			\
		 U[_idx_n],DUr[_idx_n],DUc[_idx_n],k,0,_S[_idx_n1],_S[idx],DUr,DUc); \
	}								\
	else								\
	{								\
	  solFct(_idx_n,pot[_idx_n],_un2,U[idx],			\
		 U[_idx_n],DUr[_idx_n],DUc[_idx_n],k,1,_S[_idx_n2],_S[idx],DUr,DUc); \
	}								\
      }									\
      else								\
      {								\
	  solFct(_idx_n,pot[_idx_n],U[_idx_n1],U[idx],			\
		 U[_idx_n],DUr[_idx_n],DUc[_idx_n],k,0,_S[_idx_n1],_S[idx],DUr,DUc); \
      }									\
    }									\
    else								\
    {									\
      if (_S[_idx_n2] > S_BORDER)					\
      {								\
	solFct(_idx_n,pot[_idx_n],U[_idx_n2],U[idx],			\
	       U[_idx_n],DUr[_idx_n],DUc[_idx_n],k,1,_S[_idx_n2],_S[idx],DUr,DUc); \
      }									\
      else								\
      {								        \
	_unew = U[idx] + pot[_idx_n];					\
	if (_unew < U[_idx_n])						\
	{								\
	  U[_idx_n] = _unew;						\
	  _S[_idx_n] = -_S[idx];					\
	  if (ODD(k)) { DUr[_idx_n] = pot[_idx_n] * signR[_kvar][1]; DUc[_idx_n] = 0; } \
	  else { DUr[_idx_n] = 0; DUc[_idx_n] = pot[_idx_n] * signC[_kvar][1]; } \
	  _hp->update(_idx_n);						\
	}								\
      }									\
    }
//     // ------------------------------------------------------------
#define FM_UPDATE_TRIAL_DU_P(simplexFct,solFct,signR,signC,k)		\
    _idx_n1 = simplexFct(_idx_n,_kvar,0);				\
    _idx_n2 = simplexFct(_idx_n,_kvar,1);				\
    _unew = INFPLUS;							\
    if (_S[_idx_n1] > S_BORDER)						\
    {									\
      if (_S[_idx_n2] > S_BORDER)					\
      {								        \
	_un1 = U[_idx_n1]; _un2 = U[_idx_n2];				\
	if (_un1 < _un2) solFct(_idx_n,pot[_idx_n],_un1,U[idx],		\
				U[_idx_n],DUr[_idx_n],DUc[_idx_n],k,0,_idx_n1,idx); \
	else solFct(_idx_n,pot[_idx_n],_un2,U[idx],			\
		    U[_idx_n],DUr[_idx_n],DUc[_idx_n],k,1,_idx_n2,idx);	\
      }									\
      else solFct(_idx_n,pot[_idx_n],U[_idx_n1],U[idx],			\
		  U[_idx_n],DUr[_idx_n],DUc[_idx_n],k,0,_idx_n1,idx);	\
    }									\
    else								\
    {									\
      if (_S[_idx_n2] > S_BORDER)					\
	solFct(_idx_n,pot[_idx_n],U[_idx_n2],U[idx],			\
	       U[_idx_n],DUr[_idx_n],DUc[_idx_n],k,1,_idx_n2,idx);	\
      else								\
      {								        \
	_unew = U[idx] + pot[_idx_n];					\
	if (_unew < U[_idx_n])						\
	{								\
	  U[_idx_n] = _unew;						\
	  _P1[_idx_n] = idx; _P2[_idx_n] = -1;				\
	  if (ODD(k)) { DUr[_idx_n] = pot[_idx_n] * signR[_kvar][1]; DUc[_idx_n] = 0; } \
	  else { DUr[_idx_n] = 0; DUc[_idx_n] = pot[_idx_n] * signC[_kvar][1]; } \
	  _hp->update(_idx_n);						\
	}								\
      }									\
    }
    // ------------------------------------------------------------
#define FM_UPDATE_FAR(dst,rule)			\
    rule;					\
    U[_idx_n] = U[idx] + dst;			\
    _hp->push(_idx_n);
    // ------------------------------------------------------------
#define FM_UPDATE_FAR_DU(dst,rule,signR,signC,k)			\
    FM_UPDATE_FAR(dst,rule)						\
    if (ODD(k)) { DUr[_idx_n] = pot[_idx_n] * signR[_kvar][1]; DUc[_idx_n] = 0; } \
    else { DUr[_idx_n] = 0; DUc[_idx_n] = pot[_idx_n] * signC[_kvar][1]; }
    // ------------------------------------------------------------
#define FM_UPDATE_FAR_DU_8DIAG(dst,rule)				\
    FM_UPDATE_FAR(dst,rule)						\
    DUr[_idx_n] = pot[_idx_n] * SIGNR8[_kvar][1]; DUc[_idx_n] = pot[_idx_n] * SIGNC8[_kvar][1];
    // ------------------------------------------------------------
#define FM_UPDATE(condFar,condTrial)					\
    _st = _S[_idx_n];							\
    if (_st == S_FAR)							\
    {									\
      condFar								\
      continue;								\
    }									\
    if (_st == S_TRIAL)							\
    {									\
      condTrial								\
    }
     // ------------------------------------------------------------
#define FM_UPDATE_SV(condFar,condTrial)					\
    _st = _S[_idx_n];							\
    if (_st == S_FAR)							\
    {									\
      condFar								\
      continue;								\
    }									\
    if (_st < S_FAR)							\
    {									\
      condTrial								\
    }
    // ============================================================
    // update neighbors, distance only
    // ============================================================
    template <class Grd>
    void updateNeighbors(const typename Grd::Adj4 &idx, FT *pot, FT *U)
    {
      for (_kvar = 0; _kvar < 4; _kvar++)
      {
	_idx_n = _G.neighbor4(idx,_kvar);
	FM_UPDATE(FM_UPDATE_FAR(pot[_idx_n],_S[_idx_n]=S_TRIAL),
		  FM_UPDATE_TRIAL(simplexNeighbor4,fmSolution4))
      }
    }
    // ------------------------------------------------------------
    template <class Grd>
    void updateNeighbors(const typename Grd::Adj8 &idx, FT *pot, FT *U)
    {
      FOREACH_4_8_NEIGHBOR(_kvar)
      {
	_idx_n = _G.neighbor8(idx,_kvar);
	FM_UPDATE(FM_UPDATE_FAR(pot[_idx_n],_S[_idx_n]=S_TRIAL),
		  FM_UPDATE_TRIAL(simplexNeighbor8,fmSolution8))
      }
      FOREACH_8DIAG_NEIGHBOR(_kvar)
      {
	_idx_n = _G.neighbor8(idx,_kvar);
	_st = _S[_idx_n];
	if (_st == S_FAR) { FM_UPDATE_FAR(_SQRT2*pot[_idx_n],_S[_idx_n]=S_TRIAL) }
      }
    }
    // ============================================================
    // update neighbors, distance and gradient
    // ============================================================
    template <class Grd>
    void updateNeighbors(const typename Grd::Adj4 &idx, FT *pot, FT *U, FT *DUr, FT *DUc)
    {
      for (_kvar = 0; _kvar < 4; _kvar++)
      {
	_idx_n = _G.neighbor4(idx,_kvar);
	FM_UPDATE(FM_UPDATE_FAR_DU(pot[_idx_n],_S[_idx_n]=S_TRIAL,SIGNR4,SIGNC4,_kvar),
		  FM_UPDATE_TRIAL_DU(simplexNeighbor4,fmSolution4,SIGNR4,SIGNC4,_kvar))
      }
    }
    // ------------------------------------------------------------
    template <class Grd>
    void updateNeighbors(const typename Grd::Adj8 &idx, FT *pot, FT *U, FT *DUr, FT *DUc)
    {
      IT m = 0;
      for (_kvar = 0; _kvar < 7; _kvar += 2, m++)
      {
	_idx_n = _G.neighbor8(idx,_kvar);
	FM_UPDATE(FM_UPDATE_FAR_DU(pot[_idx_n],_S[_idx_n]=S_TRIAL,SIGNR8,SIGNC8,m),
		  FM_UPDATE_TRIAL_DU(simplexNeighbor8,fmSolution8,SIGNR8,SIGNC8,m))
      }
      FOREACH_8DIAG_NEIGHBOR(_kvar)
      {
	_idx_n = _G.neighbor8(idx,_kvar);
	_st = _S[_idx_n];
	if (_st == S_FAR) { FM_UPDATE_FAR_DU_8DIAG(_SQRT2*pot[_idx_n],_S[_idx_n]=S_TRIAL) }
      }
    }
    // ============================================================
    // update neighbors, distance and gradient, and predecessor maps
    // ============================================================
    template <class Grd>
    void updateNeighborsPred(const typename Grd::Adj4 &idx, FT *pot, FT *U, FT *DUr, FT *DUc)
    {
      for (_kvar = 0; _kvar < 4; _kvar++)
      {
	_idx_n = _G.neighbor4(idx,_kvar);
	FM_UPDATE(FM_UPDATE_FAR_DU(pot[_idx_n],_S[_idx_n]=S_TRIAL;_P1[_idx_n]=idx;_P2[_idx_n]=-1,SIGNR4,SIGNC4,_kvar),
		  FM_UPDATE_TRIAL_DU_P(simplexNeighbor4,fmSolution4,SIGNR4,SIGNC4,_kvar))
      }
    }
    // ------------------------------------------------------------
    template <class Grd>
    void updateNeighborsPred(const typename Grd::Adj8 &idx, FT *pot, FT *U, FT *DUr, FT *DUc)
    {
      IT m = 0;
      for (_kvar = 0; _kvar < 7; _kvar += 2, m++)
      {
	_idx_n = _G.neighbor8(idx,_kvar);
	FM_UPDATE(FM_UPDATE_FAR_DU(pot[_idx_n],_S[_idx_n]=S_TRIAL;_P1[_idx_n]=idx;_P2[_idx_n]=-1,SIGNR8,SIGNC8,m),
		  FM_UPDATE_TRIAL_DU_P(simplexNeighbor8,fmSolution8,SIGNR8,SIGNC8,m))
      }
      FOREACH_8DIAG_NEIGHBOR(_kvar)
      {
	_idx_n = _G.neighbor8(idx,_kvar);
	_st = _S[_idx_n];
	if (_st == S_FAR) { FM_UPDATE_FAR_DU_8DIAG(_SQRT2*pot[_idx_n],_S[_idx_n]=S_TRIAL;_P1[_idx_n]=idx;_P2[_idx_n]=-1;) }
      }
    }
    // ============================================================
    // update labeled neighbors
    // ============================================================
    template <class Grd>
    bool updateLabelizedNeighbors(const typename Grd::Adj4 &idx, FT *pot, FT *U, FT *DUr, FT *DUc)
    {
      bool res = false;
      _label = _S[idx];
      for (_kvar = 0; _kvar < 4; _kvar++)
      {
	_idx_n = _G.neighbor4(idx,_kvar);
	FM_UPDATE_SV(FM_UPDATE_FAR_DU(pot[_idx_n],_S[_idx_n]=-_S[idx],SIGNR4,SIGNC4,_kvar),
		     FM_UPDATE_TRIAL_DU_SV(simplexNeighbor4,fmSolution4,SIGNR4,SIGNC4,_kvar) continue;)
      }
      return res;
    }
    // ------------------------------------------------------------
    template <class Grd, class DelGraph>
    bool updateLabelizedNeighbors(const typename Grd::Adj4 &idx, FT *pot, FT *U, FT *DUr, FT *DUc, DelGraph &dg)
    {
      bool res = false;
      _label = _S[idx];
      for (_kvar = 0; _kvar < 4; _kvar++)
      {
	_idx_n = _G.neighbor4(idx,_kvar);
	FM_UPDATE_SV(FM_UPDATE_FAR_DU(pot[_idx_n],_S[_idx_n]=-_S[idx],SIGNR4,SIGNC4,_kvar),
		     FM_UPDATE_TRIAL_DU_SV(simplexNeighbor4,fmSolution4,SIGNR4,SIGNC4,_kvar) continue;)
	// detect front meeting
	if (_st > S_BORDER && _st != _label)
	{
	  res = dg.newEdge(_st,_label,_idx_n,idx,(U[_idx_n]+U[idx])+(pot[idx]+pot[_idx_n])/2.0);
	}
      }
      FOREACH_8DIAG_NEIGHBOR(_kvar)
      {
	_idx_n = _G.neighbor8(idx,_kvar);
	_st = _S[_idx_n];
	// detect front meeting
	if (_st > S_BORDER && _st != _label)
	{
	  res = dg.newEdge(_st,_label,_idx_n,idx,(U[_idx_n]+U[idx])+(pot[idx]+pot[_idx_n])*(_SQRT2/2.0));
	}
      }
      return res;
    }
    // ------------------------------------------------------------
    template <class Grd>
    bool updateLabelizedNeighbors(const typename Grd::Adj8 &idx, FT *pot, FT *U, FT *DUr, FT *DUc)
    {
      IT m = 0;
      bool res = false;
      _label = (FT)_S[idx];
      for (_kvar = 0; _kvar < 7; _kvar += 2, m++)
      {
	_idx_n = _G.neighbor8(idx,_kvar);
	FM_UPDATE_SV(FM_UPDATE_FAR_DU(pot[_idx_n],_S[_idx_n]=-_S[idx],SIGNR8,SIGNC8,m),
		     FM_UPDATE_TRIAL_DU_SV(simplexNeighbor8,fmSolution8,SIGNR8,SIGNC8,m) continue;)
      }
      FOREACH_8DIAG_NEIGHBOR(_kvar)
      {
	_idx_n = _G.neighbor8(idx,_kvar);
	_st = _S[_idx_n];
	if (_st == S_FAR) { FM_UPDATE_FAR_DU_8DIAG(_SQRT2*pot[_idx_n],_S[_idx_n]=-_S[idx]) continue; }
      }
      return res;
    }
    // ------------------------------------------------------------
    template <class Grd, class DelGraph>
    bool updateLabelizedNeighbors(const typename Grd::Adj8 &idx, FT *pot, FT *U, FT *DUr, FT *DUc, DelGraph &dg)
    {
      IT m = 0;
      bool res = false;
      _label = _S[idx];
      for (_kvar = 0; _kvar < 7; _kvar += 2, m++)
      {
	_idx_n = _G.neighbor8(idx,_kvar);
	FM_UPDATE_SV(FM_UPDATE_FAR_DU(pot[_idx_n],_S[_idx_n]=-_S[idx],SIGNR8,SIGNC8,m),
		     FM_UPDATE_TRIAL_DU_SV(simplexNeighbor8,fmSolution8,SIGNR8,SIGNC8,m) continue;)
	// detect front meeting
	if (_st > S_BORDER && _st != _label)
	{
	  res = dg.newEdge(_st,_label,_idx_n,idx,(U[_idx_n]+U[idx])+(pot[idx]+pot[_idx_n])/2.0);
	}
      }
      FOREACH_8DIAG_NEIGHBOR(_kvar)
      {
	_idx_n = _G.neighbor8(idx,_kvar);
	_st = _S[_idx_n];
	if (_st == S_FAR) { FM_UPDATE_FAR_DU_8DIAG(_SQRT2*pot[_idx_n],_S[_idx_n]=-_S[idx]) continue; }
	// detect front meeting
	if (_st > S_BORDER && _st != _label)
	{
	  res = dg.newEdge(_st,_label,_idx_n,idx,(U[_idx_n]+U[idx])+(pot[idx]+pot[_idx_n])*(_SQRT2/2.0));
	}
      }
      return res;
    }
    // ============================================================
    // update with voronoi map
    // ============================================================
    template <class Grd>
    void updateNeighbors(const typename Grd::Adj4 &idx, FT *pot, FT *U, FT *DUr, FT *DUc, IT *V)
    {
      for (_kvar = 0; _kvar < 4; _kvar++)
      {
	_idx_n = _G.neighbor4(idx,_kvar);
	FM_UPDATE(FM_UPDATE_FAR_DU(pot[_idx_n],_S[_idx_n]=S_TRIAL;V[_idx_n]=V[idx],SIGNR4,SIGNC4,_kvar),
		  FM_UPDATE_TRIAL_DU_V(simplexNeighbor4,fmSolution4,SIGNR4,SIGNC4,_kvar,V))
      }
    }
    // ------------------------------------------------------------
    template <class Grd>
    void updateNeighbors(const typename Grd::Adj8 &idx, FT *pot, FT *U, FT *DUr, FT *DUc, IT *V)
    {
      IT m = 0;
      for (_kvar = 0; _kvar < 7; _kvar += 2, m++)
      {
	_idx_n = _G.neighbor8(idx,_kvar);
	FM_UPDATE(FM_UPDATE_FAR_DU(pot[_idx_n],_S[_idx_n]=S_TRIAL;V[_idx_n]=V[idx],SIGNR8,SIGNC8,m),
		  FM_UPDATE_TRIAL_DU_V(simplexNeighbor8,fmSolution8,SIGNR8,SIGNC8,m,V))
      }
      FOREACH_8DIAG_NEIGHBOR(_kvar)
      {
	_idx_n = _G.neighbor8(idx,_kvar);
	_st = _S[_idx_n];
	if (_st == S_FAR) { FM_UPDATE_FAR_DU_8DIAG(_SQRT2*pot[_idx_n],_S[_idx_n]=S_TRIAL;V[_idx_n]=V[idx]) }
      }
    }
    // ============================================================
    // update
    // ============================================================
    template <class Grd, class WeightClass>
    void updateNeighborsDirected(const typename Grd::Adj4 &idx, const WeightClass &W, FT *U, FT *DUr, FT *DUc)
    {
      // for (_kvar = 0; _kvar < 4; _kvar++)
      // {
      // 	_idx_n = _G.neighbor4(idx,_kvar);
      // 	FM_UPDATE(FM_UPDATE_FAR_DU(pot[_idx_n],_S[_idx_n]=S_TRIAL,SIGNR4,SIGNC4,_kvar),
      // 		  FM_UPDATE_TRIAL_DU(simplexNeighbor4,fmSolution4,SIGNR4,SIGNC4,_kvar))
      // }
    }
    // ------------------------------------------------------------
    template <class Grd, class DelGraph>
    void updateNeighbors(const typename Grd::Adj4 &idx, FT *pot, FT *U, FT *DUr, FT *DUc, IT *V, DelGraph &DG)
    {
      IT lab = V[idx], labn = -1, i1, i2, idx8;
      FT v12, v1, v2, u_idx = U[idx], u_min = INFPLUS;
      bool f = false, eorn = false;
      typename DelGraph::HEClass::HEDataClass::HEClass *he = NULL, *het = NULL;
      typename DelGraph::HEClass::HEDataClass::NodeClass *nd = NULL, *ndt = NULL;
      for (_kvar = 0; _kvar < 4; _kvar++)
      {
	idx8 = _G.neighbor8(idx,2*_kvar+1);
	_st = _S[idx8];
	if (_st > S_BORDER && V[idx8] != lab)
	{
	  labn = V[idx8];
	  i1 = idx+DOT8[_kvar][0];
	  v1 = U[idx8];//(u_idx + U[idx8]) + //(pot[idx]+pot[idx8])/_SQRT2;
	    //(pot[i1]+pot[i1+1]+pot[i1+_G.rows()]+pot[i1+_G.rows()+1])*_SQRT2/4.0;
	  typename DelGraph::HEClass *dhe = DG.isNeighbor(lab,labn);
	  if (dhe == NULL) dhe = DG.newEdge(lab,labn);
	  ndt = dhe->data->insert(i1,v1,idx,idx8);
	  //if (v1 < u_min) { u_min = v1; f = true; eorn = false; nd = ndt; }
	}
	_idx_n = _G.neighbor4(idx,_kvar);
	FM_UPDATE(FM_UPDATE_FAR_DU(pot[_idx_n],_S[_idx_n]=S_TRIAL;V[_idx_n]=V[idx],SIGNR4,SIGNC4,_kvar),
		  FM_UPDATE_TRIAL_DU_V(simplexNeighbor4,fmSolution4,SIGNR4,SIGNC4,_kvar,V);continue;)
	if (_st > S_BORDER && V[_idx_n] != lab)
	{
	  labn = V[_idx_n];
	  i1 = idx+DOT4[_kvar][0];
	  i2 = idx+DOT4[_kvar][1];
	  v12 = U[_idx_n];//(u_idx + U[_idx_n]) + (pot[_idx_n]+pot[idx])/2.0;
	  v1 = U[_idx_n];//(u_idx + U[_idx_n]) + INFPLUS;//(pot[_idx_n]+pot[idx])/_SQRT2;//(pot[i1]+pot[i1+1]+pot[i1+_G.rows()]+pot[i1+_G.rows()+1])*_SQRT2/4.0;
	  v2 = U[_idx_n];//(u_idx + U[_idx_n]) + INFPLUS;//(pot[_idx_n]+pot[idx])/_SQRT2;//(pot[i2]+pot[i2+1]+pot[i2+_G.rows()]+pot[i2+_G.rows()+1])*_SQRT2/4.0;
	  typename DelGraph::HEClass *dhe = DG.isNeighbor(lab,labn);
	  if (dhe == NULL) dhe = DG.newEdge(lab,labn);
	  het = dhe->data->insert(i1,i2,v1,v2,v12,idx,_idx_n);
	  //if (v12 < u_min) { u_min = v12; f = true; eorn = true; he = het; }
	}
      }
      /*if (f)
      {
	if (eorn) { he->data->keep = true; }
	else nd->data->keep = true;
	}*/
      // FOREACH_8DIAG_NEIGHBOR(_kvar)
      // {
      // 	_idx_n = _G.neighbor8(idx,_kvar);
      // 	_st = _S[_idx_n];
      // 	if (_st > S_BORDER && V[_idx_n] != lab)
      // 	{
      // 	  labn = V[_idx_n];
      // 	  u_tmp = (u_idx + U[_idx_n]) + (pot[_idx_n]+pot[idx])/_SQRT2;
      // 	  typename DelGraph::HEClass *he = DG.isNeighbor(lab,labn);
      // 	  if (he == NULL) DG.newEdge(lab,labn,u_tmp,idx,_idx_n);
      // 	  else
      // 	  {
      // 	    if (u_tmp < he->data->u)
      // 	    {
      // 	      he->data->u = u_tmp;
      // 	      he->data->s1 = idx;
      // 	      he->data->s2 = _idx_n;
      // 	    }
      // 	  }
      // 	}
      // }
    }
    // ------------------------------------------------------------
    template <class Grd, class WeightClass>
    void updateNeighborsDirected(const typename Grd::Adj8 &idx, const WeightClass &W, FT *U, FT *DUr, FT *DUc)
    {
      IT m = 0;
      for (_kvar = 0; _kvar < 7; _kvar += 2, m++)
      {
	_idx_n = _G.neighbor8(idx,_kvar);
	_st = _S[_idx_n];
	if (_st == S_FAR)
	{
	  U[_idx_n] = U[idx] + W(idx,_kvar);
	  _S[_idx_n] = S_TRIAL;
	  if (ODD(m)) { DUr[_idx_n] = W(idx,_kvar) * SIGNR8[_kvar][1]; DUc[_idx_n] = 0; }
	  else { DUr[_idx_n] = 0; DUc[_idx_n] = W(idx,_kvar) * SIGNC8[_kvar][1]; }
	  _hp->push(_idx_n);
	  continue;
	}
	if (_st == S_TRIAL)
	{
	  _idx_n1 = simplexNeighbor8(_idx_n,_kvar,0);
	  _idx_n2 = simplexNeighbor8(_idx_n,_kvar,1);
	  _unew = INFPLUS;
	  if (_S[_idx_n1] > S_BORDER)
	  {
	    _un1 = U[_idx_n1];
	    fmSolution8(_idx_n,W(_idx_n1,_wsn8[0][_kvar]),W(idx,_kvar),_un1,U[idx],U[_idx_n],DUr[_idx_n],DUc[_idx_n],m,0);
	  }
	  if (_S[_idx_n2] > S_BORDER)
	  {
	    _un2 = U[_idx_n2];
	    fmSolution8(_idx_n,W(_idx_n2,_wsn8[1][_kvar]),W(idx,_kvar),_un2,U[idx],U[_idx_n],DUr[_idx_n],DUc[_idx_n],m,1);
	  }
	}
      }
      FOREACH_8DIAG_NEIGHBOR(_kvar)
      {
	_idx_n = _G.neighbor8(idx,_kvar);
	_st = _S[_idx_n];
	if (_st == S_FAR)
	{
	  U[_idx_n] = U[idx] + _SQRT2*W(idx,_kvar);
	  _S[_idx_n] = S_TRIAL;
	  DUr[_idx_n] = W(idx,_kvar) * SIGNR8[_kvar][1]; DUc[_idx_n] = W(idx,_kvar) * SIGNC8[_kvar][1];
	  _hp->push(_idx_n);
	}
      }
    }
    // ------------------------------------------------------------
    template <class Grd, class DelGraph>
    void updateNeighbors(const typename Grd::Adj8 &idx, FT *pot, FT *U, FT *DUr, FT *DUc, IT *V, DelGraph &DG)
    {
      IT m = 0, lab = V[idx], labn = -1, i1, i2;
      FT v1, v2, v12, u_idx = U[idx];
      for (_kvar = 0; _kvar < 7; _kvar += 2, m++)
      {
	_idx_n = _G.neighbor8(idx,_kvar);
	FM_UPDATE(FM_UPDATE_FAR_DU(pot[_idx_n],_S[_idx_n]=S_TRIAL;V[_idx_n]=V[idx],SIGNR8,SIGNC8,m),
		  FM_UPDATE_TRIAL_DU_V(simplexNeighbor8,fmSolution8,SIGNR8,SIGNC8,m,V);continue;)
	if (_st > S_BORDER && V[_idx_n] != lab)
	{
	  labn = V[_idx_n];
	  i1 = idx+DOT4[_kvar][0];
	  i2 = idx+DOT4[_kvar][1];
	  v12 = (u_idx + U[_idx_n]) + (pot[_idx_n]+pot[idx])/2.0;
	  v1 = (u_idx + U[_idx_n]) + (pot[i1]+pot[i1+1]+pot[i1+_G.rows()]+pot[i1+_G.rows()+1])*_SQRT2/4.0;
	  v2 = (u_idx + U[_idx_n]) + (pot[i2]+pot[i2+1]+pot[i2+_G.rows()]+pot[i2+_G.rows()+1])*_SQRT2/4.0;
	  typename DelGraph::HEClass *he = DG.isNeighbor(lab,labn);
	  if (he == NULL) he = DG.newEdge(lab,labn);
	  he->data->insert(i1,i2,v1,v2,v12,idx,_idx_n);
	}
      }
      FOREACH_8DIAG_NEIGHBOR(_kvar)
      {
	_idx_n = _G.neighbor8(idx,_kvar);
	_st = _S[_idx_n];
	if (_st == S_FAR) { FM_UPDATE_FAR_DU_8DIAG(_SQRT2*pot[_idx_n],_S[_idx_n]=S_TRIAL;V[_idx_n]=V[idx]) continue; }
	// if (_st > S_BORDER && V[_idx_n] != lab)
	// {
	//   labn = V[_idx_n];
	//   i1 = idx+DOT82[_kvar][0];
	//   v1 = (u_idx + U[_idx_n]) + (pot[i1]+pot[i1+1]+pot[i1+_G.rows()]+pot[i1+_G.rows()+1])*_SQRT2/4.0;
	//   typename DelGraph::HEClass *he = DG.isNeighbor(lab,labn);
	//   if (he == NULL) he = DG.newEdge(lab,labn);
	//   he->data->insert(i1,v1,idx,_idx_n);
	// }
	// if (_st > S_BORDER && V[_idx_n] != lab)
	// {
	//   labn = V[_idx_n];
	//   u_tmp = (u_idx + U[_idx_n]) + (pot[_idx_n]+pot[idx])/_SQRT2;
	//   typename DelGraph::HEClass *he = DG.isNeighbor(lab,labn);
	//   if (he == NULL) DG.newEdge(lab,labn,u_tmp,idx,_idx_n);
	//   else
	//   {
	//     if (u_tmp < he->data->u)
	//     {
	//       he->data->u = u_tmp;
	//       he->data->s1 = idx;
	//       he->data->s2 = _idx_n;
	//     }
	//   }
	// }
      }
    }
    // ------------------------------------------------------------
    /**
     * @brief Test if a vertex is inside the current front, that is if it has no neighbor vertex S_TRIAL or S_FAR (only accepted neighbors or boundary neighbors).
     * @param[in] pos Index of the vertex to test.
     */
    template <class Grd>
    bool isInsideFront(const typename Grd::Adj4 &idx)
    {
      IT n_idx;
      unsigned short k;
      FOREACH_4_NEIGHBOR(k)
      {
	n_idx = _G.neighbor4(idx,k);
	if (_S[n_idx] < 0) return false;
      }
      return true;
    }
    // ------------------------------------------------------------
    template <class Grd>
    bool isInsideFront(const typename Grd::Adj8 &idx)
    {
      IT n_idx;
      unsigned short k;
      FOREACH_4_NEIGHBOR(k)
      {
	n_idx = _G.neighbor4(idx,k);
	if (_S[n_idx] < 0) return false;
      }
      FOREACH_8DIAG_NEIGHBOR(k)
      {
	n_idx = _G.neighbor8(idx,k);
	if (_S[n_idx] < 0) return false;
      }
      return true;
    }
    // ============================================================
    // updateNeighbors inside and outside methods
    // ============================================================
    template <class Grd>
    void updateNeighborsInside(const typename Grd::Adj4 &idx, FT *pot, FT *U)
    {
      bool t = true;
      FOREACH_4_NEIGHBOR(_kvar)
      {
	_idx_n = _G.neighbor4(idx,_kvar);
	FM_UPDATE(FM_UPDATE_FAR(pot[_idx_n],_S[_idx_n]=S_TRIAL);t=false;,
		  FM_UPDATE_TRIAL(simplexNeighbor4,fmSolution4);t=false;continue;)
	if (_st == S_DEAD_INSIDE)
	{
	  if (isInsideFront<Grd>(_idx_n))
	  {
	    _nbAcceptedInside--;
	    _S[_idx_n] = S_DEAD;
	  }
	  continue;
	}
	if ((_st == S_DEAD_OUTSIDE) && (isInsideFront<Grd>(_idx_n))) _S[_idx_n] = S_DEAD;
      }
      if (t) { _nbAcceptedInside--; _S[idx] = S_DEAD; }
    }
    // ------------------------------------------------------------
    template <class Grd>
    void updateNeighborsOutside(const typename Grd::Adj4 &idx, FT *pot, FT *U)
    {
      bool t = true;
      FOREACH_4_NEIGHBOR(_kvar)
      {
	_idx_n = _G.neighbor4(idx,_kvar);
	FM_UPDATE(FM_UPDATE_FAR(pot[_idx_n],_S[_idx_n]=S_TRIAL);t=false;,
		  FM_UPDATE_TRIAL(simplexNeighbor4,fmSolution4);t=false;continue;)
	if (_st == S_DEAD_INSIDE)
	{
	  if (isInsideFront<Grd>(_idx_n))
	  {
	    _nbAcceptedInside--;
	    _S[_idx_n] = S_DEAD;
	  }
	  continue;
	}
	if ((_st == S_DEAD_OUTSIDE) && (isInsideFront<Grd>(_idx_n))) _S[_idx_n] = S_DEAD;
      }
      if (t) { _S[idx] = S_DEAD; }
    }
    // ------------------------------------------------------------
    template <class Grd>
    void updateNeighborsInside(const typename Grd::Adj8 &idx, FT *pot, FT *U)
    {
      bool t = true;
      for (_kvar = 0; _kvar < 7; _kvar += 2)
      {
	_idx_n = _G.neighbor8(idx,_kvar);
	FM_UPDATE(FM_UPDATE_FAR(pot[_idx_n],_S[_idx_n]=S_TRIAL);t=false;,
		  FM_UPDATE_TRIAL(simplexNeighbor8,fmSolution8);t=false;continue;)
	if (_st == S_DEAD_INSIDE)
	{
	  if (isInsideFront<Grd>(_idx_n))
	  {
	    _nbAcceptedInside--;
	    _S[_idx_n] = S_DEAD;
	  }
	  continue;
	}
	if ((_st == S_DEAD_OUTSIDE) && (isInsideFront<Grd>(_idx_n))) _S[_idx_n] = S_DEAD;
      }
      FOREACH_8DIAG_NEIGHBOR(_kvar)
      {
	_idx_n = _G.neighbor8(idx,_kvar);
	_st = _S[_idx_n];
	if (_st == S_FAR) { FM_UPDATE_FAR(_SQRT2*pot[_idx_n],_S[_idx_n]=S_TRIAL); t=false; continue; }
	if (_st == S_DEAD_INSIDE)
	{
	  if (isInsideFront<Grd>(_idx_n))
	  {
	    _nbAcceptedInside--;
	    _S[_idx_n] = S_DEAD;
	  }
	  continue;
	}
	if ((_st == S_DEAD_OUTSIDE) && (isInsideFront<Grd>(_idx_n))) _S[_idx_n] = S_DEAD;
      }
      if (t) { _nbAcceptedInside--; _S[idx] = S_DEAD; }
    }
    // ------------------------------------------------------------
    template <class Grd>
    void updateNeighborsOutside(const typename Grd::Adj8 &idx, FT *pot, FT *U)
    {
      IT m = 0;
      bool t = true;
      for (_kvar = 0; _kvar < 7; _kvar += 2, m++)
      {
	_idx_n = _G.neighbor8(idx,_kvar);
	FM_UPDATE(FM_UPDATE_FAR(pot[_idx_n],_S[_idx_n]=S_TRIAL); t=false;,
		  FM_UPDATE_TRIAL(simplexNeighbor8,fmSolution8); t=false;continue;)
	if (_st == S_DEAD_INSIDE)
	{
	  if (isInsideFront<Grd>(_idx_n))
	  {
	    _nbAcceptedInside--;
	    _S[_idx_n] = S_DEAD;
	  }
	  continue;
	}
	if ((_st == S_DEAD_OUTSIDE) && (isInsideFront<Grd>(_idx_n))) _S[_idx_n] = S_DEAD;
      }
      FOREACH_8DIAG_NEIGHBOR(_kvar)
      {
	_idx_n = _G.neighbor8(idx,_kvar);
	_st = _S[_idx_n];
	if (_st == S_FAR) { FM_UPDATE_FAR(_SQRT2*pot[_idx_n],_S[_idx_n]=S_TRIAL); t=false; continue; }
	if (_st == S_DEAD_INSIDE)
	{
	  if (isInsideFront<Grd>(_idx_n))
	  {
	    _nbAcceptedInside--;
	    _S[_idx_n] = S_DEAD;
	  }
	  continue;
	}
	if ((_st == S_DEAD_OUTSIDE) && (isInsideFront<Grd>(_idx_n))) _S[_idx_n] = S_DEAD;
      }
      if (t) { _S[idx] = S_DEAD; }
    }
    // ============================================================
    // updateNeighborsInside with gradient
    // ============================================================
    template <class Grd>
    void updateNeighborsInside(const typename Grd::Adj4 &idx, FT *pot, FT *U, FT *DUr, FT *DUc)
    {
      bool t = true;
      FOREACH_4_NEIGHBOR(_kvar)
      {
	_idx_n = _G.neighbor4(idx,_kvar);
	FM_UPDATE(FM_UPDATE_FAR_DU(pot[_idx_n],_S[_idx_n]=S_TRIAL,SIGNR4,SIGNC4,_kvar);t=false;,
		  FM_UPDATE_TRIAL_DU(simplexNeighbor4,fmSolution4,SIGNR4,SIGNC4,_kvar);t=false;continue;)
	if (_st == S_DEAD_INSIDE)
	{
	  if (isInsideFront<Grd>(_idx_n))
	  {
	    _nbAcceptedInside--;
	    _S[_idx_n] = S_DEAD;
	  }
	  continue;
	}
	if ((_st == S_DEAD_OUTSIDE) && (isInsideFront<Grd>(_idx_n))) _S[_idx_n] = S_DEAD;
      }
      if (t) { _nbAcceptedInside--; _S[idx] = S_DEAD; }
    }
    // ------------------------------------------------------------
    template <class Grd>
    void updateNeighborsOutside(const typename Grd::Adj4 &idx, FT *pot, FT *U, FT *DUr, FT *DUc)
    {
      bool t = true;
      FOREACH_4_NEIGHBOR(_kvar)
      {
	_idx_n = _G.neighbor4(idx,_kvar);
	FM_UPDATE(FM_UPDATE_FAR_DU(pot[_idx_n],_S[_idx_n]=S_TRIAL,SIGNR4,SIGNC4,_kvar);t=false;,
		  FM_UPDATE_TRIAL_DU(simplexNeighbor4,fmSolution4,SIGNR4,SIGNC4,_kvar);t=false;continue;)
	if (_st == S_DEAD_INSIDE)
	{
	  if (isInsideFront<Grd>(_idx_n))
	  {
	    _nbAcceptedInside--;
	    _S[_idx_n] = S_DEAD;
	  }
	  continue;
	}
	if ((_st == S_DEAD_OUTSIDE) && (isInsideFront<Grd>(_idx_n))) _S[_idx_n] = S_DEAD;
      }
      if (t) { _S[idx] = S_DEAD; }
    }
    // ------------------------------------------------------------
    template <class Grd>
    void updateNeighborsInside(const typename Grd::Adj8 &idx, FT *pot, FT *U, FT *DUr, FT *DUc)
    {
      bool t = true;
      IT m = 0;
      for (_kvar = 0; _kvar < 7; _kvar += 2, m++)
      {
	_idx_n = _G.neighbor8(idx,_kvar);
	FM_UPDATE(FM_UPDATE_FAR_DU(pot[_idx_n],_S[_idx_n]=S_TRIAL,SIGNR8,SIGNC8,m);t=false;,
		  FM_UPDATE_TRIAL_DU(simplexNeighbor8,fmSolution8,SIGNR8,SIGNC8,m);t=false;continue;)
	if (_st == S_DEAD_INSIDE)
	{
	  if (isInsideFront<Grd>(_idx_n))
	  {
	    _nbAcceptedInside--;
	    _S[_idx_n] = S_DEAD;
	  }
	  continue;
	}
	if ((_st == S_DEAD_OUTSIDE) && (isInsideFront<Grd>(_idx_n))) _S[_idx_n] = S_DEAD;
      }
      FOREACH_8DIAG_NEIGHBOR(_kvar)
      {
	_idx_n = _G.neighbor8(idx,_kvar);
	_st = _S[_idx_n];
	if (_st == S_FAR) { FM_UPDATE_FAR_DU_8DIAG(_SQRT2*pot[_idx_n],_S[_idx_n]=S_TRIAL) t=false; continue; }
	if (_st == S_DEAD_INSIDE)
	{
	  if (isInsideFront<Grd>(_idx_n))
	  {
	    _nbAcceptedInside--;
	    _S[_idx_n] = S_DEAD;
	  }
	  continue;
	}
	if ((_st == S_DEAD_OUTSIDE) && (isInsideFront<Grd>(_idx_n))) _S[_idx_n] = S_DEAD;
      }
      if (t) { _nbAcceptedInside--; _S[idx] = S_DEAD; }
    }
    // ------------------------------------------------------------
    template <class Grd>
    void updateNeighborsOutside(const typename Grd::Adj8 &idx, FT *pot, FT *U, FT *DUr, FT *DUc)
    {
      bool t = true;
      IT m = 0;
      for (_kvar = 0; _kvar < 7; _kvar += 2, m++)
      {
	_idx_n = _G.neighbor8(idx,_kvar);
	FM_UPDATE(FM_UPDATE_FAR_DU(pot[_idx_n],_S[_idx_n]=S_TRIAL,SIGNR8,SIGNC8,m);t=false;,
		  ;		  FM_UPDATE_TRIAL_DU(simplexNeighbor8,fmSolution8,SIGNR8,SIGNC8,m);t=false;continue;)
	if (_st == S_DEAD_INSIDE)
	{
	  if (isInsideFront<Grd>(_idx_n))
	  {
	    _nbAcceptedInside--;
	    _S[_idx_n] = S_DEAD;
	  }
	  continue;
	}
	if ((_st == S_DEAD_OUTSIDE) && (isInsideFront<Grd>(_idx_n))) _S[_idx_n] = S_DEAD;
      }
      FOREACH_8DIAG_NEIGHBOR(_kvar)
      {
	_idx_n = _G.neighbor8(idx,_kvar);
	_st = _S[_idx_n];
	if (_st == S_FAR) { FM_UPDATE_FAR_DU_8DIAG(_SQRT2*pot[_idx_n],_S[_idx_n]=S_TRIAL);t=false; continue; }
	if (_st == S_DEAD_INSIDE)
	{
	  if (isInsideFront<Grd>(_idx_n))
	  {
	    _nbAcceptedInside--;
	    _S[_idx_n] = S_DEAD;
	  }
	  continue;
	}
	if ((_st == S_DEAD_OUTSIDE) && (isInsideFront<Grd>(_idx_n))) _S[_idx_n] = S_DEAD;
      }
      if (t) { _S[idx] = S_DEAD; }
    }
    // ============================================================
    // updateNeighborsInside with gradient and Voronoi maps
    // ============================================================
    template <class Grd>
    void updateNeighborsInside(const typename Grd::Adj4 &idx, FT *pot, FT *U, FT *DUr, FT *DUc, IT *V, const IT &cur_lab)
    {
      bool t = true;
      FOREACH_4_NEIGHBOR(_kvar)
      {
	_idx_n = _G.neighbor4(idx,_kvar);
	FM_UPDATE(FM_UPDATE_FAR_DU(pot[_idx_n],_S[_idx_n]=S_TRIAL,SIGNR4,SIGNC4,_kvar);t=false;,
		  FM_UPDATE_TRIAL_DU(simplexNeighbor4,fmSolution4,SIGNR4,SIGNC4,_kvar);t=false;continue;)
	if (_st == S_DEAD_INSIDE)
	{
	  if (isInsideFront<Grd>(_idx_n))
	  {
	    _nbAcceptedInside--;
	    _S[_idx_n] = S_DEAD;
	    V[_idx_n] = cur_lab;
	  }
	  continue;
	}
	if ((_st == S_DEAD_OUTSIDE) && (isInsideFront<Grd>(_idx_n))) _S[_idx_n] = S_DEAD;
      }
      if (t) { _nbAcceptedInside--; _S[idx] = S_DEAD; V[idx] = cur_lab; }
    }
    // ------------------------------------------------------------
    template <class Grd>
    void updateNeighborsOutside(const typename Grd::Adj4 &idx, FT *pot, FT *U, FT *DUr, FT *DUc, IT *V, const IT &cur_lab)
    {
      bool t = true;
      FOREACH_4_NEIGHBOR(_kvar)
      {
	_idx_n = _G.neighbor4(idx,_kvar);
	FM_UPDATE(FM_UPDATE_FAR_DU(pot[_idx_n],_S[_idx_n]=S_TRIAL,SIGNR4,SIGNC4,_kvar);t=false;,
		  FM_UPDATE_TRIAL_DU(simplexNeighbor4,fmSolution4,SIGNR4,SIGNC4,_kvar);t=false;continue;)
	if (_st == S_DEAD_INSIDE)
	{
	  if (isInsideFront<Grd>(_idx_n))
	  {
	    _nbAcceptedInside--;
	    _S[_idx_n] = S_DEAD;
	    V[_idx_n] = cur_lab;
	  }
	  continue;
	}
	if ((_st == S_DEAD_OUTSIDE) && (isInsideFront<Grd>(_idx_n))) _S[_idx_n] = S_DEAD;
      }
      if (t) { _S[idx] = S_DEAD; }
    }
    // ------------------------------------------------------------
    template <class Grd>
    void updateNeighborsInside(const typename Grd::Adj8 &idx, FT *pot, FT *U, FT *DUr, FT *DUc, IT *V, const IT &cur_lab)
    {
      bool t = true;
      IT m = 0;
      for (_kvar = 0; _kvar < 7; _kvar += 2, m++)
      {
	_idx_n = _G.neighbor8(idx,_kvar);
	FM_UPDATE(FM_UPDATE_FAR_DU(pot[_idx_n],_S[_idx_n]=S_TRIAL,SIGNR8,SIGNC8,m);t=false;,
		  FM_UPDATE_TRIAL_DU(simplexNeighbor8,fmSolution8,SIGNR8,SIGNC8,m);t=false;continue;)
	if (_st == S_DEAD_INSIDE)
	{
	  if (isInsideFront<Grd>(_idx_n))
	  {
	    _nbAcceptedInside--;
	    _S[_idx_n] = S_DEAD;
	    V[_idx_n] = cur_lab;
	  }
	  continue;
	}
	if ((_st == S_DEAD_OUTSIDE) && (isInsideFront<Grd>(_idx_n))) _S[_idx_n] = S_DEAD;
      }
      FOREACH_8DIAG_NEIGHBOR(_kvar)
      {
	_idx_n = _G.neighbor8(idx,_kvar);
	_st = _S[_idx_n];
	if (_st == S_FAR) { FM_UPDATE_FAR_DU_8DIAG(_SQRT2*pot[_idx_n],_S[_idx_n]=S_TRIAL) t=false; continue; }
	if (_st == S_DEAD_INSIDE)
	{
	  if (isInsideFront<Grd>(_idx_n))
	  {
	    _nbAcceptedInside--;
	    _S[_idx_n] = S_DEAD;
	    V[_idx_n] = cur_lab;
	  }
	  continue;
	}
	if ((_st == S_DEAD_OUTSIDE) && (isInsideFront<Grd>(_idx_n))) _S[_idx_n] = S_DEAD;
      }
      if (t) { _nbAcceptedInside--; _S[idx] = S_DEAD; V[idx] = cur_lab; }
    }
    // ------------------------------------------------------------
    template <class Grd>
    void updateNeighborsOutside(const typename Grd::Adj8 &idx, FT *pot, FT *U, FT *DUr, FT *DUc, IT *V, const IT &cur_lab)
    {
      IT m = 0;
      bool t = true;
      for (_kvar = 0; _kvar < 7; _kvar += 2, m++)
      {
	_idx_n = _G.neighbor8(idx,_kvar);
	FM_UPDATE(FM_UPDATE_FAR_DU(pot[_idx_n],_S[_idx_n]=S_TRIAL,SIGNR8,SIGNC8,m);t=false;,
		  FM_UPDATE_TRIAL_DU(simplexNeighbor8,fmSolution8,SIGNR8,SIGNC8,m);t=false;continue;)
	if (_st == S_DEAD_INSIDE)
	{
	  if (isInsideFront<Grd>(_idx_n))
	  {
	    _nbAcceptedInside--;
	    _S[_idx_n] = S_DEAD;
	    V[_idx_n] = cur_lab;
	  }
	  continue;
	}
	if ((_st == S_DEAD_OUTSIDE) && (isInsideFront<Grd>(_idx_n))) _S[_idx_n] = S_DEAD;
      }
      FOREACH_8DIAG_NEIGHBOR(_kvar)
      {
	_idx_n = _G.neighbor8(idx,_kvar);
	_st = _S[_idx_n];
	if (_st == S_FAR) { FM_UPDATE_FAR_DU_8DIAG(_SQRT2*pot[_idx_n],_S[_idx_n]=S_TRIAL);t=false; continue; }
	if (_st == S_DEAD_INSIDE)
	{
	  if (isInsideFront<Grd>(_idx_n))
	  {
	    _nbAcceptedInside--;
	    _S[_idx_n] = S_DEAD;
	    V[_idx_n] = cur_lab;
	  }
	  continue;
	}
	if ((_st == S_DEAD_OUTSIDE) && (isInsideFront<Grd>(_idx_n))) _S[_idx_n] = S_DEAD;
      }
      if (t) { _S[idx] = S_DEAD; }
    }
    // =======================================================================
    // updateNeighborsInside with gradient and Voronoi maps and Delaunay graph
    // =======================================================================
    template <class Grd, class DelGraph>
    bool updateNeighborsInside(const typename Grd::Adj4 &idx, FT *pot, FT *U, FT *DUr, FT *DUc, IT *V, const IT &cur_lab, DelGraph &DG, FT *Upred)
    {
      bool t = true, res = false;
      FOREACH_4_NEIGHBOR(_kvar)
      {
	_idx_n = _G.neighbor4(idx,_kvar);
	FM_UPDATE(FM_UPDATE_FAR_DU(pot[_idx_n],_S[_idx_n]=S_TRIAL,SIGNR4,SIGNC4,_kvar);t=false;,
		  FM_UPDATE_TRIAL_DU(simplexNeighbor4,fmSolution4,SIGNR4,SIGNC4,_kvar);t=false;continue;)
	if (_st == S_DEAD_INSIDE)
	{
	  if (isInsideFront<Grd>(_idx_n))
	  {
	    _nbAcceptedInside--;
	    _S[_idx_n] = S_DEAD;
	    V[_idx_n] = cur_lab;
	  }
	  continue;
	}
	if (_st == S_DEAD_OUTSIDE)
	{
	  if (isInsideFront<Grd>(_idx_n)) _S[_idx_n] = S_DEAD;
	  res = DG.newEdge(V[_idx_n],cur_lab,_idx_n,idx,(U[_idx_n]+U[idx])+(pot[idx]+pot[_idx_n])/2.0);
	}
      }
      if (t) { _nbAcceptedInside--; _S[idx] = S_DEAD; V[idx] = cur_lab; }
      return res;
    }
    // ------------------------------------------------------------
    template <class Grd, class DelGraph>
    bool updateNeighborsOutside(const typename Grd::Adj4 &idx, FT *pot, FT *U, FT *DUr, FT *DUc, IT *V, const IT &cur_lab, DelGraph &DG, FT *Upred)
    {
      bool t = true, res = false;
      FOREACH_4_NEIGHBOR(_kvar)
      {
	_idx_n = _G.neighbor4(idx,_kvar);
	FM_UPDATE(FM_UPDATE_FAR_DU(pot[_idx_n],_S[_idx_n]=S_TRIAL,SIGNR4,SIGNC4,_kvar);t=false;,
		  FM_UPDATE_TRIAL_DU(simplexNeighbor4,fmSolution4,SIGNR4,SIGNC4,_kvar);t=false;continue;)
	if (_st == S_DEAD_INSIDE)
	{
	  if (isInsideFront<Grd>(_idx_n))
	  {
	    _nbAcceptedInside--;
	    _S[_idx_n] = S_DEAD;
	    V[_idx_n] = cur_lab;
	  }
	  continue;
	}
	if ((_st == S_DEAD_OUTSIDE) && (isInsideFront<Grd>(_idx_n))) _S[_idx_n] = S_DEAD;
      }
      if (t) { _S[idx] = S_DEAD; }
      return res;
    }
    // ------------------------------------------------------------
    template <class Grd, class DelGraph>
    bool updateNeighborsInside(const typename Grd::Adj8 &idx, FT *pot, FT *U, FT *DUr, FT *DUc, IT *V, const IT &cur_lab, DelGraph &DG, FT *Upred)
    {
      bool t = true, res = false;
      IT m = 0;
      for (_kvar = 0; _kvar < 7; _kvar += 2, m++)
      {
	_idx_n = _G.neighbor8(idx,_kvar);
	FM_UPDATE(FM_UPDATE_FAR_DU(pot[_idx_n],_S[_idx_n]=S_TRIAL,SIGNR8,SIGNC8,m);t=false;,
		  FM_UPDATE_TRIAL_DU(simplexNeighbor8,fmSolution8,SIGNR8,SIGNC8,m);t=false;continue;)
	if (_st == S_DEAD_INSIDE)
	{
	  if (isInsideFront<Grd>(_idx_n))
	  {
	    _nbAcceptedInside--;
	    _S[_idx_n] = S_DEAD;
	    V[_idx_n] = cur_lab;
	  }
	  continue;
	}
	if (_st == S_DEAD_OUTSIDE)
	{
	  if (isInsideFront<Grd>(_idx_n)) _S[_idx_n] = S_DEAD;
	  res = DG.newEdge(V[_idx_n],cur_lab,_idx_n,idx,(Upred[_idx_n]+U[idx])+(pot[idx]+pot[_idx_n])/2.0);
	  continue;
	}
	if (_st == S_DEAD && V[_idx_n] != cur_lab)
	  res = DG.newEdge(V[_idx_n],cur_lab,_idx_n,idx,(Upred[_idx_n]+U[idx])+(pot[idx]+pot[_idx_n])/2.0);
      }
      FOREACH_8DIAG_NEIGHBOR(_kvar)
      {
	_idx_n = _G.neighbor8(idx,_kvar);
	_st = _S[_idx_n];
	if (_st == S_FAR) { FM_UPDATE_FAR_DU_8DIAG(_SQRT2*pot[_idx_n],_S[_idx_n]=S_TRIAL) t=false; continue; }
	if (_st == S_DEAD_INSIDE)
	{
	  if (isInsideFront<Grd>(_idx_n))
	  {
	    _nbAcceptedInside--;
	    _S[_idx_n] = S_DEAD;
	    V[_idx_n] = cur_lab;
	  }
	  continue;
	}
	if (_st == S_DEAD_OUTSIDE)
	{
	  if (isInsideFront<Grd>(_idx_n)) _S[_idx_n] = S_DEAD;
	  res = DG.newEdge(V[_idx_n],cur_lab,_idx_n,idx,(Upred[_idx_n]+U[idx])+(pot[idx]+pot[_idx_n])*(_SQRT2/2.0));
	  continue;
	}
	if (_st == S_DEAD && V[_idx_n] != cur_lab)
	  res = DG.newEdge(V[_idx_n],cur_lab,_idx_n,idx,(Upred[_idx_n]+U[idx])+(pot[idx]+pot[_idx_n])*(_SQRT2/2.0));
      }
      if (t) { _nbAcceptedInside--; _S[idx] = S_DEAD; V[idx] = cur_lab; }
      return res;
    }
    // ------------------------------------------------------------
    template <class Grd, class DelGraph>
    bool updateNeighborsOutside(const typename Grd::Adj8 &idx, FT *pot, FT *U, FT *DUr, FT *DUc, IT *V, const IT &cur_lab, DelGraph &DG, FT *Upred)
    {
      IT m = 0;
      bool t = true, res = false;
      for (_kvar = 0; _kvar < 7; _kvar += 2, m++)
      {
	_idx_n = _G.neighbor8(idx,_kvar);
	FM_UPDATE(FM_UPDATE_FAR_DU(pot[_idx_n],_S[_idx_n]=S_TRIAL,SIGNR8,SIGNC8,m);t=false;,
		  FM_UPDATE_TRIAL_DU(simplexNeighbor8,fmSolution8,SIGNR8,SIGNC8,m);t=false;continue;)
	if (_st == S_DEAD_INSIDE)
	{
	  if (isInsideFront<Grd>(_idx_n))
	  {
	    _nbAcceptedInside--;
	    _S[_idx_n] = S_DEAD;
	    V[_idx_n] = cur_lab;
	  }
	  res = DG.newEdge(cur_lab,V[idx],_idx_n,idx,(U[_idx_n]+Upred[idx])+(pot[idx]+pot[_idx_n])/2.0);
	  continue;
	}
	if ((_st == S_DEAD_OUTSIDE) && (isInsideFront<Grd>(_idx_n))) { _S[_idx_n] = S_DEAD; continue; }
	if (_st == S_DEAD && V[_idx_n] == cur_lab)
	  res = DG.newEdge(cur_lab,V[idx],_idx_n,idx,(U[_idx_n]+Upred[idx])+(pot[idx]+pot[_idx_n])/2.0);
      }
      FOREACH_8DIAG_NEIGHBOR(_kvar)
      {
	_idx_n = _G.neighbor8(idx,_kvar);
	_st = _S[_idx_n];
	if (_st == S_FAR) { FM_UPDATE_FAR_DU_8DIAG(_SQRT2*pot[_idx_n],_S[_idx_n]=S_TRIAL);t=false; continue; }
	if (_st == S_DEAD_INSIDE)
	{
	  if (isInsideFront<Grd>(_idx_n))
	  {
	    _nbAcceptedInside--;
	    _S[_idx_n] = S_DEAD;
	    V[_idx_n] = cur_lab;
	  }
	  res = DG.newEdge(cur_lab,V[idx],_idx_n,idx,(U[_idx_n]+Upred[idx])+(pot[idx]+pot[_idx_n])*(_SQRT2/2.0));
	  continue;
	}
	if ((_st == S_DEAD_OUTSIDE) && (isInsideFront<Grd>(_idx_n))) { _S[_idx_n] = S_DEAD; continue; }
	if (_st == S_DEAD && V[_idx_n] == cur_lab)
	  res = DG.newEdge(cur_lab,V[idx],_idx_n,idx,(U[_idx_n]+Upred[idx])+(pot[idx]+pot[_idx_n])*(_SQRT2/2.0));
      }
      if (t) { _S[idx] = S_DEAD; }
      return res;
    }
    // ============================================================
    //
    // ============================================================
    // template <class Grd>
    // void updateNeighborsInside(const typename Grd::Adj4 &idx, FT *pot, FT *U, FT *DUr, FT *DUc, IT *V)
    // {
    //   bool t = true, b = false;
    //   IT idx8, idx1, idx2, idxtt;
    //   FT ndut, ndu = INFPLUS;

    //   FOREACH_4_NEIGHBOR(_kvar)
    //   {
    // 	_idx_n = _G.neighbor4(idx,_kvar);

    // 	if (_S[_idx_n] == S_BORDER) { b = true; continue; }

    // 	idx8 = _G.neighbor8(idx,2*_kvar+1);
    // 	if ((_S[idx8] == S_DEAD_OUTSIDE) || (_S[idx8] == S_DEAD && V[_idx_n] != _label))
    // 	{
    // 	  ndut = U[idx8];
    // 	  if (_S[_idx_n] == S_DEAD_OUTSIDE || (_S[_idx_n] == S_DEAD && V[_idx_n] != _label))
    // 	  {
    // 	    idx1 = _medialgraph->addNode(idx+DOT82[_kvar][0]);
    // 	    idx2 = _medialgraph->addNode(idx+DOT82[_kvar][1]);
    // 	    if (_medialgraph->isNeighbor(idx1,idx2) < 0) _medialgraph->link(idx1,idx2);
    // 	  }
    // 	  idxtt = _G.neighbor4(idx,(_kvar+1)%4);
    // 	  if (_S[idxtt] == S_DEAD_OUTSIDE || (_S[idxtt] == S_DEAD && V[idxtt] != _label))
    // 	  {
    // 	    idx1 = _medialgraph->addNode(idx+DOT82[_kvar][0]);
    // 	    idx2 = _medialgraph->addNode(idx+DOT82[_kvar][2]);
    // 	    if (_medialgraph->isNeighbor(idx1,idx2) < 0) _medialgraph->link(idx1,idx2);
    // 	  }
    // 	}

    // 	_idx_n = _G.neighbor4(idx,_kvar);
    // 	FM_UPDATE(FM_UPDATE_FAR_DU(pot[_idx_n],_S[_idx_n]=S_TRIAL,SIGNR4,SIGNC4,_kvar);t=false;_P1[_idx_n]=idx;_P2[_idx_n]=-1;,
    // 		  FM_UPDATE_TRIAL_DU_P(simplexNeighbor4,fmSolution4,SIGNR4,SIGNC4,_kvar);t=false;continue;)
    // 	if (_st == S_DEAD_INSIDE)
    // 	{
    // 	  if (isInsideFront<Grd>(_idx_n))
    // 	  {
    // 	    _nbAcceptedInside--;
    // 	    _S[_idx_n] = S_DEAD;
    // 	  }
    // 	  continue;
    // 	}
    // 	if ((_st == S_DEAD_OUTSIDE) && (isInsideFront<Grd>(_idx_n))) _S[_idx_n] = S_DEAD;
    //   }
    //   if (t) { _nbAcceptedInside--; _S[idx] = S_DEAD; }
    // }
    // ------------------------------------------------------------
    template <class Grd>
    void updateNeighborsInside(const typename Grd::Adj8 &idx, FT *pot, FT *U, FT *DUr, FT *DUc, IT *V)
    {
      bool t = true;
      IT m = 0;
      for (_kvar = 0; _kvar < 7; _kvar += 2, m++)
      {
	_idx_n = _G.neighbor8(idx,_kvar);
	FM_UPDATE(FM_UPDATE_FAR_DU(pot[_idx_n],_S[_idx_n]=S_TRIAL,SIGNR8,SIGNC8,m);t=false;,
		  FM_UPDATE_TRIAL_DU(simplexNeighbor8,fmSolution8,SIGNR8,SIGNC8,m);t=false;continue;)
	if (_st == S_DEAD_INSIDE)
	{
	  if (isInsideFront<Grd>(_idx_n))
	  {
	    _nbAcceptedInside--;
	    _S[_idx_n] = S_DEAD;
	  }
	  continue;
	}
	if ((_st == S_DEAD_OUTSIDE) && (isInsideFront<Grd>(_idx_n))) _S[_idx_n] = S_DEAD;
      }
      FOREACH_8DIAG_NEIGHBOR(_kvar)
      {
	_idx_n = _G.neighbor8(idx,_kvar);
	_st = _S[_idx_n];
	if (_st == S_FAR) { FM_UPDATE_FAR_DU_8DIAG(_SQRT2*pot[_idx_n],_S[_idx_n]=S_TRIAL) t=false; continue; }
	if (_st == S_DEAD_INSIDE)
	{
	  if (isInsideFront<Grd>(_idx_n))
	  {
	    _nbAcceptedInside--;
	    _S[_idx_n] = S_DEAD;
	  }
	  continue;
	}
	if ((_st == S_DEAD_OUTSIDE) && (isInsideFront<Grd>(_idx_n))) _S[_idx_n] = S_DEAD;
      }
      if (t) _nbAcceptedInside--;
    }
    // ------------------------------------------------------------
    // template <class Grd>
    // void updateNeighborsOutside(const typename Grd::Adj4 &idx, FT *pot, FT *U, FT *DUr, FT *DUc, FT *Upred, FT *DUr_prev, FT *DUc_prev, IT *V, FT *DUS)
    // {
    //   FT ndu = INFPLUS, ndut = INFPLUS;
    //   bool t = true, f = false, b = false, dorn = false;
    //   IT idxk = -1, idx8, kp, idx1, idx2, idxtt, dt;

    //   FOREACH_4_NEIGHBOR(_kvar)
    //   {
    // 	_idx_n = _G.neighbor4(idx,_kvar);

    // 	if (_S[_idx_n] == S_BORDER) { b = true; continue; }

    // 	idx8 = _G.neighbor8(idx,2*_kvar+1);
    // 	if (V[idx8] == _label)
    // 	{
    // 	  idx1 = -1;
    // 	  ndut = U[idx8] + Upred[idx8];
    // 	  if (_S[_idx_n] == S_DEAD_OUTSIDE || (_S[_idx_n] == S_DEAD && V[_idx_n] != _label))
    // 	  {
    // 	    idx1 = _medialgraph->addNode(idx+DOT8[_kvar][0],(U[idx]+Upred[idx] + ndut)/4.0,idx,idx8);
    // 	    idx2 = _medialgraph->addNode(idx+DOT8[_kvar][1]);
    // 	    if (_medialgraph->isNeighbor(idx1,idx2) < 0) _medialgraph->link(idx1,idx2);
    // 	  }
    // 	  idxtt = _G.neighbor4(idx,(_kvar+1)%4);
    // 	  if (_S[idxtt] == S_DEAD_OUTSIDE || (_S[idxtt] == S_DEAD && V[idxtt] != _label))
    // 	  {
    // 	    idx1 = _medialgraph->addNode(idx+DOT8[_kvar][0],(U[idx]+Upred[idx] + ndut)/4.0,idx,idx8);
    // 	    idx2 = _medialgraph->addNode(idx+DOT8[_kvar][2]);
    // 	    if (_medialgraph->isNeighbor(idx1,idx2) < 0) _medialgraph->link(idx1,idx2);
    // 	  }
    // 	  if (ndut < ndu && idx1 > 0) { ndu = ndut; idxk = idx1; dorn = true; f = true; }
    // 	}

    // 	FM_UPDATE(FM_UPDATE_FAR_DU(pot[_idx_n],_S[_idx_n]=S_TRIAL,SIGNR4,SIGNC4,_kvar);t=false;_P1[_idx_n]=idx;_P2[_idx_n]=-1;,
    // 		  FM_UPDATE_TRIAL_DU_P(simplexNeighbor4,fmSolution4,SIGNR4,SIGNC4,_kvar);t=false;continue;)

    // 	if (_st == S_DEAD_INSIDE)
    // 	{
    // 	  if (isInsideFront<Grd>(_idx_n))
    // 	  {
    // 	    _nbAcceptedInside--;
    // 	    _S[_idx_n] = S_DEAD;
    // 	  }
    // 	  ndut = U[_idx_n] + Upred[_idx_n];
    // 	  idx1 = _medialgraph->addNode(idx+DOT4[_kvar][0]);
    // 	  idx2 = _medialgraph->addNode(idx+DOT4[_kvar][1]);
    // 	  if ((dt = _medialgraph->isNeighbor(idx1,idx2)) < 0) dt = _medialgraph->link(idx1,idx2,(U[idx]+Upred[idx] + ndut)/4.0,idx,_idx_n);
    // 	  else _medialgraph->updateDart(dt,(U[idx]+Upred[idx] + ndut)/4.0,idx,_idx_n);
    // 	  if (ndut < ndu) { ndu = ndut; idxk = dt; f = true; dorn = false; }
    // 	  continue;
    // 	}
    // 	if (_st == S_DEAD_OUTSIDE)
    // 	{
    // 	  if (isInsideFront<Grd>(_idx_n)) _S[_idx_n] = S_DEAD;
    // 	  continue;
    // 	}
    // 	if (_st == S_DEAD)
    // 	{
    // 	  if (V[_idx_n] == _label)
    // 	  {
    // 	    ndut = U[_idx_n] + Upred[_idx_n];
    // 	    idx1 = _medialgraph->addNode(idx+DOT4[_kvar][0]);
    // 	    idx2 = _medialgraph->addNode(idx+DOT4[_kvar][1]);
    // 	    if ((dt = _medialgraph->isNeighbor(idx1,idx2)) < 0) dt = _medialgraph->link(idx1,idx2,(U[idx]+Upred[idx] + ndut)/4.0,idx,_idx_n);
    // 	    else _medialgraph->updateDart(dt,(U[idx]+Upred[idx] + ndut)/4.0,idx,_idx_n);
    // 	    if (ndut < ndu) { ndu = ndut; idxk = dt; f = true; dorn = false; }
    // 	  }
    // 	}
    //   }
    //   if (t) { _S[idx] = S_DEAD; }

    //   if (f) // point of the medial set
    //   {
    // 	if (dorn) _medialgraph->node(idxk)->keep = true;
    // 	else { _medialgraph->dart(idxk)->keep = true; _medialgraph->dart(_medialgraph->dart(idxk)->alpha)->keep = true; }
    // 	DUS[idx] = (U[idx] + Upred[idx] + ndu) / 4.0;
    //   }
    // }
    // ------------------------------------------------------------
    template <class Grd>
    void updateNeighborsOutside(const typename Grd::Adj8 &idx, FT *pot, FT *U, FT *DUr, FT *DUc, FT *Uprev, FT *DUr_prev, FT *DUc_prev, IT *V, FT *DUS)
    {
      IT m = 0;
      FT n, nprev, dur, duc, durp, ducp, ndu, ndut = INFPLUS, nduk = INFPLUS, dus = -1;
      typename Grd::Adj8 idxk = -1;

      for (_kvar = 0; _kvar < 7; _kvar += 2, m++)
      {
	_idx_n = _G.neighbor8(idx,_kvar);
	FM_UPDATE(FM_UPDATE_FAR_DU(pot[_idx_n],_S[_idx_n]=S_TRIAL,SIGNR8,SIGNC8,m),
		  FM_UPDATE_TRIAL_DU(simplexNeighbor8,fmSolution8,SIGNR8,SIGNC8,m);continue;)
	if (_st == S_DEAD_INSIDE)
	{
	  if (isInsideFront<Grd>(_idx_n))
	  {
	    _nbAcceptedInside--;
	    _S[_idx_n] = S_DEAD;
	  }

	  // nprev = std::sqrt(DUr_prev[_idx_n]*DUr_prev[_idx_n] + DUc_prev[_idx_n]*DUc_prev[_idx_n]);
	  // durp = DUr_prev[_idx_n]/nprev;
	  // ducp = DUc_prev[_idx_n]/nprev;

	  // n = std::sqrt(DUr[_idx_n]*DUr[_idx_n] + DUc[_idx_n]*DUc[_idx_n]);
	  // dur = DUr[_idx_n]/n;
	  // duc = DUc[_idx_n]/n;

	  // ndut = (std::abs(dur+durp) + std::abs(duc+ducp))/2.;
	  // if (ndut < DUS[_idx_n]) DUS[_idx_n] = ndut;
	  // idxk = 1;

	  continue;
	}
	if ((_st == S_DEAD_OUTSIDE) && (isInsideFront<Grd>(_idx_n))) _S[_idx_n] = S_DEAD;
      }
      FOREACH_8DIAG_NEIGHBOR(_kvar)
      {
	_idx_n = _G.neighbor8(idx,_kvar);
	_st = _S[_idx_n];
	if (_st == S_FAR) { FM_UPDATE_FAR_DU_8DIAG(_SQRT2*pot[_idx_n],_S[_idx_n]=S_TRIAL); continue; }
	if (_st == S_DEAD_INSIDE)
	{
	  if (isInsideFront<Grd>(_idx_n))
	  {
	    _nbAcceptedInside--;
	    _S[_idx_n] = S_DEAD;
	  }

	  // nprev = std::sqrt(DUr_prev[_idx_n]*DUr_prev[_idx_n] + DUc_prev[_idx_n]*DUc_prev[_idx_n]);
	  // durp = DUr_prev[_idx_n]/nprev;
	  // ducp = DUc_prev[_idx_n]/nprev;

	  // n = std::sqrt(DUr[_idx_n]*DUr[_idx_n] + DUc[_idx_n]*DUc[_idx_n]);
	  // dur = DUr[_idx_n]/n;
	  // duc = DUc[_idx_n]/n;

	  // ndut = (std::abs(dur+durp) + std::abs(duc+ducp))/2.;
	  // if (ndut < DUS[_idx_n]) DUS[_idx_n] = ndut;
	  // idxk = 1;
	  continue;
	}
	if ((_st == S_DEAD_OUTSIDE) && (isInsideFront<Grd>(_idx_n))) _S[_idx_n] = S_DEAD;
      }

      // if (idxk >= 0) // case of a medialset pair (idx,idxk) with minimal gradient sum
      // {
	// nprev = std::sqrt(DUr_prev[idx]*DUr_prev[idx] + DUc_prev[idx]*DUc_prev[idx]);
	// durp = DUr_prev[idx]/nprev;
	// ducp = DUc_prev[idx]/nprev;

	// n = std::sqrt(DUr[idx]*DUr[idx] + DUc[idx]*DUc[idx]);
	// dur = DUr[idx]/n;
	// duc = DUc[idx]/n;

	// ndut = (std::abs(dur+durp) + std::abs(duc+ducp))/2.;
	// if (ndut < DUS[idx]) DUS[idx] = ndut;

      // }

    }
    // ------------------------------------------------------------
    /*void fmUpdateNeighbors4(const IT &idx, FT *pot, FT *U, FT *DUr, FT *DUc)
    {
      IT st;
      IT idx_n, idx_n1, idx_n2;
      FT unew, un1, un2;
      // -- update neighbors of idx -----------
      for (unsigned short k = 0; k < 4; k++)
      {
	idx_n = neighbor4(idx,k);
	st = _S[idx_n];
	// -- not yet treated ------------------
	if (st == S_FAR)
	{
	  _S[idx_n] = S_TRIAL;
	  U[idx_n] = U[idx] + pot[idx_n];
	  if (ODD(k)) { DUr[idx_n] = pot[idx_n] * SIGNR4[k][1]; DUc[idx_n] = 0; }
	  else { DUr[idx_n] = 0; DUc[idx_n] = pot[idx_n] * SIGNC4[k][1]; }
	  _hp->push(idx_n);
	  continue;
	}
	// --  already treated (in the heap), try to reduce --
	if (st == S_TRIAL)
	{
	  idx_n1 = simplexNeighbor4(idx_n,k,0);
	  idx_n2 = simplexNeighbor4(idx_n,k,1);
	  unew = INFPLUS;
	  if (_S[idx_n1] == S_DEAD)
	  {
	    if (_S[idx_n2] == S_DEAD)
	    {
	      un1 = U[idx_n1];
	      un2 = U[idx_n2];
	      if (un1 < un2) fmUpdateNeighbor4(idx_n,pot[idx_n],un1,U[idx],
					       U[idx_n],DUr[idx_n],DUc[idx_n],k,0);
	      else fmUpdateNeighbor4(idx_n,pot[idx_n],un2,U[idx],
				     U[idx_n],DUr[idx_n],DUc[idx_n],k,1);
	    }
	    else fmUpdateNeighbor4(idx_n,pot[idx_n],U[idx_n1],U[idx],
				   U[idx_n],DUr[idx_n],DUc[idx_n],k,0);
	  }
	  else // idx_n1 is not S_DEAD
	  {
	    if (_S[idx_n2] == S_DEAD)
	      fmUpdateNeighbor4(idx_n,pot[idx_n],U[idx_n2],U[idx],
				U[idx_n],DUr[idx_n],DUc[idx_n],k,1);
	    else
	    {
	      unew = U[idx] + pot[idx_n];
	      if (unew < U[idx_n])
	      {
	      U[idx_n] = unew;
		if (ODD(k)) { DUr[idx_n] = pot[idx_n] * SIGNR4[k][1]; DUc[idx_n] = 0; }
		else { DUr[idx_n] = 0; DUc[idx_n] = pot[idx_n] * SIGNC4[k][1]; }
		_hp->update(idx_n);
	      }
	    }
	  }
	}
      }
      }*/
    // ------------------------------------------------------------
    /*void fmUpdateNeighbors8(const IT &idx, FT *pot, FT *U, FT *DUr, FT *DUc)
    {
      IT st;
      IT idx_n, idx_n1, idx_n2;
      FT unew, un1, un2;
      unsigned short k = 0, m = 0;
      // -- update 4-neighbors of idx -----------
      for (; k < 7; k+=2, m++)
      {
	idx_n = neighbor8(idx,k);
	st = _S[idx_n];
	// -- not yet treated ------------------
	if (st == S_FAR)
	{
	  _S[idx_n] = S_TRIAL;
	  U[idx_n] = U[idx] + pot[idx_n];
	  if (ODD(m)) { DUr[idx_n] = pot[idx_n] * SIGNR8[k][1]; DUc[idx_n] = 0; }
	  else { DUr[idx_n] = 0; DUc[idx_n] = pot[idx_n] * SIGNC8[k][1]; }
	  _hp->push(idx_n);
	  continue;
	}
	// --  already treated (in the heap), try to reduce --
	if (st == S_TRIAL)
	{
	  idx_n1 = simplexNeighbor8(idx_n,k,0);
	  idx_n2 = simplexNeighbor8(idx_n,k,1);
	  unew = INFPLUS;
	  if (_S[idx_n1] == S_DEAD)
	  {
	    if (_S[idx_n2] == S_DEAD)
	    {
	      un1 = U[idx_n1];
	      un2 = U[idx_n2];
	      if (un1 < un2) fmUpdateNeighbor8(idx_n,pot[idx_n],un1,U[idx],
					       U[idx_n],DUr[idx_n],DUc[idx_n],k,0,m);
	      else fmUpdateNeighbor8(idx_n,pot[idx_n],un2,U[idx],
				     U[idx_n],DUr[idx_n],DUc[idx_n],k,1,m);
	    }
	    else fmUpdateNeighbor8(idx_n,pot[idx_n],U[idx_n1],U[idx],
				   U[idx_n],DUr[idx_n],DUc[idx_n],k,0,m);
	  }
	  else // idx_n1 is not S_DEAD
	  {
	    if (_S[idx_n2] == S_DEAD)
	      fmUpdateNeighbor8(idx_n,pot[idx_n],U[idx_n2],U[idx],
				U[idx_n],DUr[idx_n],DUc[idx_n],k,1,m);
	    else
	    {
	      unew = U[idx] + pot[idx_n];
	      if (unew < U[idx_n])
	      {
		U[idx_n] = unew;
		if (ODD(m)) { DUr[idx_n] = pot[idx_n] * SIGNR8[k][1]; DUc[idx_n] = 0; }
		else { DUr[idx_n] = 0; DUc[idx_n] = pot[idx_n] * SIGNC8[k][1]; }
		_hp->update(idx_n);
	      }
	    }
	  }
	}
      }
      // --- update diagonal neighbors of idx ----------------
      for (k = 1; k < 8; k+=2)
      {
	idx_n = neighbor8(idx,k);
	st = _S[idx_n];
	// -- not yet treated ------------------
	if (st == S_FAR)
	{
	  _S[idx_n] = S_TRIAL;
	  U[idx_n] = U[idx] + _SQRT2 * pot[idx_n];
	  DUr[idx_n] = pot[idx_n] * SIGNR8[k][1]; DUc[idx_n] = pot[idx_n] * SIGNC8[k][1];
	  _hp->push(idx_n);
	}
      }
      }*/
    // ------------------------------------------------------------
    /* void fmUpdateLabelizedNeighbors(const IT &idx, FT *pot, FT *U,
				    const unsigned short &adjacency,
				    const unsigned short &radius)
    {
      IT st, snew, st1, st2;
      IT idx_n, idx_n1, idx_n2;
      FT unew, un1, un2, u_idx;
      // -- update neighbors of idx -----------------
      for (unsigned short k = 0; k < adjacency; k++)
      {
	idx_n = neighbor(idx,radius,k);
	st = _S[idx_n];
	// -- not yet treated ------------------
	if (st == S_FAR)
	{
	  _S[idx_n] = -_S[idx];
	  U[idx_n] = U[idx] + _W[radius][k] * pot[idx_n];
	  _hp->push(idx_n);
	  continue;
	}
	// --  already treated (in the heap), try to reduce --
	if (st < S_FAR)
	{
	  idx_n1 = simplexNeighbor(idx_n,radius,k,0);
	  idx_n2 = simplexNeighbor(idx_n,radius,k,1);
	  unew = INFPLUS;
	  u_idx = U[idx];
	  st1 = _S[idx_n1]; st2 = _S[idx_n2];
	  if (st1 > S_BORDER)
	  {
	    if (st2 > S_BORDER)
	    {
	      un1 = U[idx_n1];
	      un2 = U[idx_n2];
	      if (un1 < un2) // from (idx,idx_1)
	      {
		unew = fmSolution(pot[idx_n],un1,u_idx,radius,k);
		snew = (ODD(k) ?
			(un1 < u_idx + radius*pot[idx_n]*(_SQRT2-1.0) ? -st1 : -_S[idx]) :
			(un1+radius*pot[idx_n]*(_SQRT2-1.0) < u_idx ? -st1 : -_S[idx]));
	      }
	      else // from (idx,idx_n2)
	      {
		unew = fmSolution(pot[idx_n],un2,u_idx,radius,k);
		snew = (ODD(k) ?
			(un2 < u_idx + radius*pot[idx_n]*(_SQRT2-1.0) ? -st2 : -_S[idx]) :
			(un2+radius*pot[idx_n]*(_SQRT2-1.0) < u_idx ? -st2 : -_S[idx]));
	      }
	    }
	    else // from (idx,idx_1)
	    {
	      unew = fmSolution(pot[idx_n],U[idx_n1],u_idx,radius,k);
	      snew = (ODD(k) ?
		      (un1 < u_idx + radius*pot[idx_n]*(_SQRT2-1.0) ? -st1 : -_S[idx]) :
		      (un1+radius*pot[idx_n]*(_SQRT2-1.0) < u_idx ? -st1 : -_S[idx]));
	    }
	  }
	  else // idx_n1 is not S_DEAD
	  {
	    if (_S[idx_n2] > S_BORDER)
	    {
	      un2 = U[idx_n2];
	      unew = fmSolution(pot[idx_n],un2,u_idx,radius,k);
	      snew = (ODD(k) ?
		      (un2 < u_idx + radius*pot[idx_n]*(_SQRT2-1.0) ? -st2 : -_S[idx]) :
		      (un2+radius*pot[idx_n]*(_SQRT2-1.0) < u_idx ? -st2 : -_S[idx]));
	    }
	    else
	    {
	      unew = u_idx + _W[radius][k] * pot[idx_n];
	      snew = -_S[idx];
	    }
	  }
	  if (unew < U[idx_n])
	  {
	    _S[idx_n] = snew;
	    U[idx_n] = unew;
	    _hp->update(idx_n);
	  }
	}
      }
      }*/
    // ============================================================
    // Propagation methods
    // ============================================================
    /**
     * @brief Propagation according to a potential map.
     * @param[in] pot Potential map.
     * @param[out] U Distance map.
     */
    void propagate(FT *pot, FT *U)
    {
      IT idx;
      while (!_hp->empty())
      {
	idx = _hp->pop();
	_S[idx] = S_DEAD;
	updateNeighbors<GridClass>(idx,pot,U);
      }
    }
    // ------------------------------------------------------------
    /**
     * @brief Propagation according to a potential map.
     * @param[in] pot Potential map.
     * @param[out] U Distance map.
     * @param[out] DUr Row gradient of the distance.
     * @param[out] DUc Column gradient of the distance.
     */
    void propagate(FT *pot, FT *U, FT *DUr, FT *DUc)
    {
      IT idx;
      while (!_hp->empty())
      {
	idx = _hp->pop();
	_S[idx] = S_DEAD;
	updateNeighbors<GridClass>(idx,pot,U,DUr,DUc);
      }
    }
    // ------------------------------------------------------------
    /**
     * @brief Propagation according to a potential map.
     * @param[in] pot Potential map.
     * @param[out] U Distance map.
     * @param[out] DUr Row gradient of the distance.
     * @param[out] DUc Column gradient of the distance.
     * @param[out] V Voronoi map.
     */
    void propagate(FT *pot, FT *U, FT *DUr, FT *DUc, IT *V)
    {
      IT idx;
      while (!_hp->empty())
      {
	idx = _hp->pop();
	_S[idx] = S_DEAD;
	updateNeighbors<GridClass>(idx,pot,U,DUr,DUc,V);
      }
    }
    // ------------------------------------------------------------
    /**
     * @brief Propagation according to a potential map.
     * @param[in] pot Potential map.
     * @param[out] U Distance map.
     * @param[out] DUr Row gradient of the distance.
     * @param[out] DUc Column gradient of the distance.
     * @param[out] V Voronoi map.
     */
    void propagatePred(FT *pot, FT *U, FT *DUr, FT *DUc)
    {
      IT idx;
      while (!_hp->empty())
      {
	idx = _hp->pop();
	_S[idx] = S_DEAD;
	updateNeighborsPred<GridClass>(idx,pot,U,DUr,DUc);
      }
    }
    // ============================================================
    // Propagation methods with Voronoi (labelized _S)
    // ============================================================
    /**
     * @brief Propagation according to a potential map.
     * @param[in] pot Potential map.
     * @param[out] U Distance map.
     * @param[out] DUr Row gradient of the distance.
     * @param[out] DUc Column gradient of the distance.
     * @param[out] V Voronoi map.
     * @param[out] DG Delaunay graph.
     */
    template <class DelGraph>
    void propagateLabeled(FT *pot, FT *U, FT *DUr, FT *DUc, DelGraph &DG)
    {
      IT idx;
      while (!_hp->empty())
      {
	idx = _hp->pop();
	_S[idx] = -_S[idx];
	updateLabelizedNeighbors<GridClass,DelGraph>(idx,pot,U,DUr,DUc,DG);
      }
    }
    // ------------------------------------------------------------
    /**
     * @brief Propagation according to a potential map.
     * @param[in] pot Potential map.
     * @param[out] U Distance map.
     * @param[out] DUr Row gradient of the distance.
     * @param[out] DUc Column gradient of the distance.
     * @param[out] V Voronoi map.
     */
    void propagateLabeled(FT *pot, FT *U, FT *DUr, FT *DUc)
    {
      IT idx;
      while (!_hp->empty())
      {
	idx = _hp->pop();
	_S[idx] = -_S[idx];
	updateLabelizedNeighbors<GridClass>(idx,pot,U,DUr,DUc);
      }
    }
    // ------------------------------------------------------------
    /**
     * @brief Propagation according to a potential map.
     * @param[in] pot Potential map.
     * @param[out] U Distance map.
     * @param[out] DUr Row gradient of the distance.
     * @param[out] DUc Column gradient of the distance.
     * @param[out] V Voronoi map.
     * @param[out] DG Delaunay graph.
     */
    template <class DelGraph>
    void propagateLabeled(FT *pot, FT *U, FT *DUr, FT *DUc, DelGraph &DG, IT m)
    {
      IT idx;
      while (!_hp->empty())
      {
	idx = _hp->pop();
	_S[idx] = -_S[idx];
	if (updateLabelizedNeighbors<GridClass,DelGraph>(idx,pot,U,DUr,DUc,DG)) m--;
	if (m == 0) break;
      }
    }
    // ============================================================
    // Propagation methods : iterative cases
    // ============================================================
    void propagate(FT *pot, FT *Uprev, FT *U)
    {
      IT idx;
      _nbAcceptedInside = 0;
      while (!_hp->empty())
      {
	idx = _hp->pop();
	if (Uprev[idx] < U[idx]) // not inside
	{
	  _S[idx] = S_DEAD_OUTSIDE;
	  updateNeighborsOutside<GridClass>(idx,pot,U);
	}
	else  // inside
	{
	  _S[idx] = S_DEAD_INSIDE;
	  _nbAcceptedInside++;
	  updateNeighborsInside<GridClass>(idx,pot,U);
	  Uprev[idx] = U[idx];
	}
	if (_nbAcceptedInside <= 0) break;  // all accepted inside are computed
      }
    }
    // ------------------------------------------------------------
    void propagate(FT *pot, FT *Uprev, FT *DUrp, FT *DUcp, FT *U, FT *DUr, FT *DUc)
    {
      IT idx;
      _nbAcceptedInside = 0;
      while (!_hp->empty())
      {
	idx = _hp->pop();
	if (Uprev[idx] < U[idx]) // not inside
	{
	  _S[idx] = S_DEAD_OUTSIDE;
	  updateNeighborsOutside<GridClass>(idx,pot,U,DUr,DUc);
	}
	else  // inside
	{
	  _S[idx] = S_DEAD_INSIDE;
	  _nbAcceptedInside++;
	  updateNeighborsInside<GridClass>(idx,pot,U,DUr,DUc);
	  Uprev[idx] = U[idx]; DUrp[idx] = DUr[idx]; DUcp[idx] = DUc[idx];
	}
	if (_nbAcceptedInside <= 0) break;  // all accepted inside are computed
      }
    }
    // ------------------------------------------------------------
    void propagate(FT *pot, IT *V, const IT &lab_idx, FT *Uprev, FT *DUrp, FT *DUcp, FT *U, FT *DUr, FT *DUc)
    {
      IT idx;
      _nbAcceptedInside = 0;
      while (!_hp->empty())
      {
	idx = _hp->pop();
	if (Uprev[idx] < U[idx]) // not inside
	{
	  _S[idx] = S_DEAD_OUTSIDE;
	  updateNeighborsOutside<GridClass>(idx,pot,U,DUr,DUc,V,lab_idx);
	}
	else  // inside
	{
	  _S[idx] = S_DEAD_INSIDE;
	  _nbAcceptedInside++;
	  updateNeighborsInside<GridClass>(idx,pot,U,DUr,DUc,V,lab_idx);
	  Uprev[idx] = U[idx]; DUrp[idx] = DUr[idx]; DUcp[idx] = DUc[idx];
	}
	if (_nbAcceptedInside <= 0) break;  // all accepted inside are computed
      }
    }
    // ------------------------------------------------------------
    template <class DelGraph>
    void propagate(FT *pot, IT *V, const IT &lab_idx, FT *Uprev, FT *DUrp, FT *DUcp, FT *U, FT *DUr, FT *DUc, DelGraph &DG)
    {
      IT idx;
      _nbAcceptedInside = 0;
      while (!_hp->empty())
      {
	idx = _hp->pop();
	if (Uprev[idx] < U[idx]) // not inside
	{
	  _S[idx] = S_DEAD_OUTSIDE;
	  updateNeighborsOutside<GridClass,DelGraph>(idx,pot,U,DUr,DUc,V,lab_idx,DG,Uprev);
	}
	else  // inside
	{
	  _S[idx] = S_DEAD_INSIDE;
	  _nbAcceptedInside++;
	  updateNeighborsInside<GridClass,DelGraph>(idx,pot,U,DUr,DUc,V,lab_idx,DG,Uprev);
	  Uprev[idx] = U[idx]; DUrp[idx] = DUr[idx]; DUcp[idx] = DUc[idx];
	}
	if (_nbAcceptedInside <= 0) break;  // all accepted inside are computed
      }
    }
    // ============================================================
    // Propagation methods
    // ============================================================
    /**
     * @brief Propagation according to a potential map.
     * @param[in] pot Potential map.
     * @param[out] U Distance map.
     * @param[out] DUr Row gradient of the distance.
     * @param[out] DUc Column gradient of the distance.
     * @param[out] V Voronoi map.
     * @param[out] DG Delaunay graph.
     */
    template <class DelGraph>
    void propagate(FT *pot, FT *U, FT *DUr, FT *DUc, IT *V, DelGraph &DG)
    {
      IT idx;
      while (!_hp->empty())
      {
	idx = _hp->pop();
	_S[idx] = S_DEAD;
	updateNeighbors<GridClass>(idx,pot,U,DUr,DUc,V,DG);
      }
    }
    // ------------------------------------------------------------
    /**
     * @brief Propagation according to a potential map.
     * @param[in] pot Potential map.
     * @param[out] U Distance map.
     * @param[out] DUr Row gradient of the distance.
     * @param[out] DUc Column gradient of the distance.
     */
    template <class WeightClass>
    void directedPropagation(WeightClass &W, FT *U, FT *DUr, FT *DUc)
    {
      std::cerr << "propagation\n";
      IT idx, nb = 0;
      while (!_hp->empty())
      {
	idx = _hp->pop();
	_S[idx] = S_DEAD;
	updateNeighborsDirected<GridClass,WeightClass>(idx,W,U,DUr,DUc);
	nb++;
      }
      std::cerr << "end propagation\n";

    }
    // ------------------------------------------------------------
    FT propagate(const IT &idx_to, FT *pot, FT *U, FT *DUr, FT *DUc)
    {
      IT idx;
      FT res = INFPLUS;
      while (!_hp->empty())
      {
	idx = _hp->pop();
	_S[idx] = S_DEAD;
	if (idx == idx_to) { res = U[idx]; break; }
	updateNeighbors<GridClass>(idx,pot,U,DUr,DUc);
      }
      while (!_hp->empty()) { idx = _hp->pop(); U[idx] = INFPLUS; }
      return res;
    }
    // ============================================================
    // Initialization methods
    // ============================================================
    void init(FT *U, IT *S)
    {
      IT pos;
      if (_S) delete[] _S;
      if (S) _S = S;
      else _S = new IT[_G.size()];
      _G.fillBoundary(_S, S_BORDER);
      for_grid_interior(_G,pos,IT)
      {
      	U[pos] = INFPLUS;
	_S[pos] = S_FAR;
      }
      if (_hp) delete _hp;
      _hp = new MinHeap<FT,IT>(_G.size(), U);
    }
    // ------------------------------------------------------------
    void init(FT *pot, FT *U)
    {
      IT pos;
      for_grid_interior(_G,pos,IT)
      {
      	U[pos] = INFPLUS;
      	_S[pos] = (pot[pos] < 0 ? S_BORDER : S_FAR);
      }
      _hp->clear();
    }
    // ------------------------------------------------------------
    /**
     * @bried Initialize arrays for propagation (distance and state maps). The potential map also allows to represent the propagation domain: pot[idx]<0 implies that idx is outside.
     * @param[in] pot Potential map.
     * @param[out] U Distance (Minimal action map) initialized to INFPLUS everywhere.
     * @param[out] S State map initialized to S_FAR for nodes inside the domain, and to S_BORDER for nodes outside.
     */
    void init(FT *pot, FT *U, IT *S)
    {
      IT pos;
      if (_S) delete[] _S;
      if (S) _S = S;
      else _S = new IT[_G.size()];
      _G.fillBoundary(_S, S_BORDER);
      for_grid_interior(_G,pos,IT)
      {
      	U[pos] = INFPLUS;
      	_S[pos] = (pot[pos] < 0 ? S_BORDER : S_FAR);
      }
      // create min heap
      if (_hp) delete _hp;
      _hp = new MinHeap<FT,IT>(_G.size(), U);
    }
    // ------------------------------------------------------------
    void init(FT *pot, FT *U, FT *DUr, FT *DUc)
    {
      IT pos;
      for_grid_interior(_G,pos,IT)
      {
      	DUr[pos] = DUc[pos] = 0.0;
	U[pos] = INFPLUS;
      	_S[pos] = (pot[pos] < 0 ? S_BORDER : S_FAR);
      }
      _hp->clear();
    }
    // ------------------------------------------------------------
    /**
     * @bried Initialize arrays for propagation (distance and state maps). The potential map also allows to represent the propagation domain: pot[idx]<0 implies that idx is outside.
     * @param[in] pot Potential map.
     * @param[out] U Distance (Minimal action map) initialized to INFPLUS everywhere.
     * @param[out] S State map initialized to S_FAR for nodes inside the domain, and to S_BORDER for nodes outside.
     */
    void init(FT *pot, FT *U, FT *DUr, FT *DUc, IT *S)
    {
      IT pos;
      if (_S) delete[] _S;
      if (S) _S = S;
      else _S = new IT[_G.size()];
      _G.fillBoundary(_S, S_BORDER);
      _G.fillBoundary(DUr, (FT)0.0);
      _G.fillBoundary(DUc, (FT)0.0);
      for_grid_interior(_G,pos,IT)
      {
      	DUr[pos] = DUc[pos] = 0.0;
	U[pos] = INFPLUS;
      	_S[pos] = (pot[pos] < 0 ? S_BORDER : S_FAR);
      }
      // create min heap
      if (_hp) delete _hp;
      _hp = new MinHeap<FT,IT>(_G.size(), U);
    }
    // ------------------------------------------------------------
    /**
     * @bried Initialize arrays for propagation (distance and state maps). The potential map also allows to represent the propagation domain: pot[idx]<0 implies that idx is outside.
     * @param[in] pot Potential map.
     * @param[out] U Distance (Minimal action map) initialized to INFPLUS everywhere.
     * @param[out] S State map initialized to S_FAR for nodes inside the domain, and to S_BORDER for nodes outside.
     */
    void initPred(FT *pot, FT *U, FT *DUr, FT *DUc, IT *S)
    {
      IT pos;
      if (_S) delete[] _S;
      if (S) _S = S;
      else _S = new IT[_G.size()];
      _G.fillBoundary(_S, S_BORDER);
      _G.fillBoundary(DUr, (FT)0.0);
      _G.fillBoundary(DUc, (FT)0.0);
      _G.fillBoundary(_P1,-1);
      _G.fillBoundary(_P2,-1);
      for_grid_interior(_G,pos,IT)
      {
      	DUr[pos] = DUc[pos] = 0.0;
	_P1[pos] = _P2[pos] = -1;
	U[pos] = INFPLUS;
      	_S[pos] = (pot[pos] < 0 ? S_BORDER : S_FAR);
      }
      // create min heap
      if (_hp) delete _hp;
      _hp = new MinHeap<FT,IT>(_G.size(), U);
    }
    // ------------------------------------------------------------
    /**
     * @bried Initialize arrays for propagation (distance and state maps). The potential map also allows to represent the propagation domain: pot[idx]<0 implies that idx is outside.
     * @param[in] pot Potential map.
     * @param[out] U Distance (Minimal action map) initialized to INFPLUS everywhere.
     * @param[out] S State map initialized to S_FAR for nodes inside the domain, and to S_BORDER for nodes outside.
     */
    void initLabelized(FT *pot, FT *U, FT *DUr, FT *DUc)
    {
      IT pos;
      _G.fillBoundary(_S, S_BORDER);
      _G.fillBoundary(DUr, (FT)0.0);
      _G.fillBoundary(DUc, (FT)0.0);
      for_grid_interior(_G,pos,IT)
      {
      	DUr[pos] = DUc[pos] = U[pos] = INFPLUS;
      	_S[pos] = (pot[pos] < 0 ? S_BORDER : S_FAR);
      }
      // create min heap
      if (_hp) delete _hp;
      _hp = new MinHeap<FT,IT>(_G.size(), U);
    }
    // ------------------------------------------------------------
    void initFPS(FT *pot, FT *U, FT *DUr, FT *DUc, FT *Upred, IT &fp)
    {
      IT pos;
      FT mx = -1.0;
      for_grid_interior(_G,pos,IT)
      {
	if (Upred[pos] > mx) { mx = Upred[pos]; fp = pos; }
      	DUr[pos] = DUc[pos] = 0.0;
	U[pos] = INFPLUS;
      	_S[pos] = (pot[pos] < 0 ? S_BORDER : S_FAR);
      }
      _hp->clear();
    }
    // ------------------------------------------------------------
    /**
     * @bried Initialize arrays for propagation (distance and state maps). The potential map also allows to represent the propagation domain: pot[idx]<0 implies that idx is outside.
     * @param[in] pot Potential map.
     * @param[out] U Distance (Minimal action map) initialized to INFPLUS everywhere.
     * @param[out] S State map initialized to S_FAR for nodes inside the domain, and to S_BORDER for nodes outside.
     */
    void init(FT *pot, FT *U, FT *DUr, FT *DUc, IT *S, IT *V, const IT &lab)
    {
      IT pos;
      if (_S) delete[] _S;
      if (S) _S = S;
      else _S = new IT[_G.size()];
      _G.fillBoundary(_S, S_BORDER);
      _G.fillBoundary(DUr, 0.0);
      _G.fillBoundary(DUc, 0.0);
      _G.fillBoundary(V,0);
      for_grid_interior(_G,pos,IT)
      {
      	DUr[pos] = DUc[pos] = 0.0;
	U[pos] = INFPLUS;
      	if (pot[pos] < 0) { _S[pos] = S_BORDER; V[pos] = 0; }
	else { _S[pos] = S_FAR; V[pos] = lab; }
      }
      // create min heap
      if (_hp) delete _hp;
      _hp = new MinHeap<FT,IT>(_G.size(), U);
    }
    // ------------------------------------------------------------
    void init(FT *pot, FT *U, IT *Pred, IT *S)
    {
      IT pos;
      if (_S) delete[] _S;
      if (S) _S = S;
      else _S = new IT[_G.size()];
      _G.fillBoundary(_S, S_BORDER);
      for_grid_interior(_G,pos,IT)
      {
	U[pos] = INFPLUS; Pred[pos] = pos;
	_S[pos] = (pot[pos] < 0 ? S_BORDER : S_FAR);
      }
      if (_hp) delete _hp;
      _hp = new MinHeap<FT,IT>(_G.size(), U);
    }
    // ------------------------------------------------------------
    void init(FT *U, IT *Pred, IT *S)
    {
      IT pos;
      if (_S) delete[] _S;
      if (S) _S = S;
      else _S = new IT[_G.size()];
      _G.fillBoundary(_S, S_BORDER);
      for_grid_interior(_G,pos,IT)
      {
	U[pos] = INFPLUS;
	Pred[pos] = pos;
	_S[pos] = S_FAR;
      }
      if (_hp) delete _hp;
      _hp = new MinHeap<FT,IT>(_G.size(), U);
    }
    // ============================================================
    // Adding nodes
    // ============================================================
    inline void addNode(const IT &idx, FT *U)
    { U[idx] = 0; _S[idx] = S_TRIAL; _hp->push(idx); }
    // ------------------------------------------------------------
    inline void addNode(const IT &idx, FT *U, FT *DUr, FT *DUc)
    {
        DUr[idx] = DUc[idx] = U[idx] = 0;

        _S[idx] = S_TRIAL;

        _hp->push(idx);
    }
    // ------------------------------------------------------------
    inline void addNodePred(const IT &idx, FT *U, FT *DUr, FT *DUc)
    { DUr[idx] = DUc[idx] = U[idx] = 0; _S[idx] = S_TRIAL; _hp->push(idx); _P2[idx] = _P1[idx] = idx; }
    // ------------------------------------------------------------
    inline void addNode(const IT &idx, const IT &label, FT *U)
    { U[idx] = 0; _S[idx] = label; _hp->push(idx); }
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
    void addNodes(IT *seeds_idx, IT nb_seeds, FT *U, FT *DUr, FT *DUc)
    {
      IT idx;
      for (IT r = 0; r < nb_seeds; r++)
      {
	idx = seeds_idx[r];
	DUr[idx] = DUc[idx] = U[idx] = 0;
	_S[idx] = S_TRIAL;
	_hp->push(idx);
      }
    }
    // ============================================================
    // Adding labeled nodes
    // ============================================================
    /**
     * @brief Add labelized grid nodes as boundary condition of the eikonal equation.
     * @param[in] seeds_idx Array of seeds represented by their indicies.
     * @param[in] nb_seeds Number of seeds in the array.
     * @param[in] U Distance (Minimal action) map.
     */
    void addLabeledNodes(IT *seeds_idx, const IT &nb_seeds, FT *U)
    {
      IT idx;
      for (IT r = 0; r < nb_seeds; r++)
      {
	idx = seeds_idx[r];
	U[idx] = 0;
	_S[idx] = S_TRIAL;
	_hp->push(idx);
	S_TRIAL--;
      }
    }
    // ------------------------------------------------------------
    /**
     * @brief Add labelized grid nodes as boundary condition of the eikonal equation.
     * @param[in] seeds_idx Array of seeds represented by their indicies.
     * @param[in] nb_seeds Number of seeds in the array.
     * @param[in] U Distance (Minimal action) map.
     */
    void addLabeledNodes(IT *seeds_idx, const IT &nb_seeds, FT *U, FT *DUr, FT *DUc)
    {
      IT idx;
      for (IT r = 0; r < nb_seeds; r++)
      {
	idx = seeds_idx[r];
	U[idx] = 0;
	DUr[idx] = 0;
	DUc[idx] = 0;
	_S[idx] = S_TRIAL;
	_hp->push(idx);
	S_TRIAL--;
      }
    }
    // ------------------------------------------------------------
    /**
     * @brief Add labelized grid nodes as boundary condition of the eikonal equation.
     * @param[in] seeds_idx Array of seeds represented by their indicies.
     * @param[in] nb_seeds Number of seeds in the array.
     * @param[in] U Distance (Minimal action) map.
     */
    void addLabeledNodes(const IT &idx1, const IT &idx2,FT *U)
    {
      U[idx1] = 0;
      _S[idx1] = S_TRIAL;
      _hp->push(idx1);
      U[idx2] = 0;
      _S[idx2] = S_TRIAL;
      _hp->push(idx2);
      S_TRIAL--;
    }
    // ============================================================
    // Adding labeled nodes
    // ============================================================
    /**
     * @brief Add labelized grid nodes as boundary condition of the eikonal equation.
     * @param[in] seeds_idx Array of seeds represented by their indicies.
     * @param[in] nb_seeds Number of seeds in the array.
     * @param[in] U Distance (Minimal action) map.
     */
    void addLabeledNodes(IT *seeds_idx, const IT &nb_seeds, FT *U, IT *V, IT label1 = 1)
    {
      IT idx;
      for (IT r = 0; r < nb_seeds; r++, label1++)
      {
	idx = seeds_idx[r];
	U[idx] = 0;
	_S[idx] = S_TRIAL;
	V[idx] = label1;
	_hp->push(idx);
      }
    }
    // ------------------------------------------------------------
    /**
     * @brief Add labelized grid nodes as boundary condition of the eikonal equation.
     * @param[in] seeds_idx Array of seeds represented by their indicies.
     * @param[in] nb_seeds Number of seeds in the array.
     * @param[in] U Distance (Minimal action) map.
     */
    void addLabeledNodes(IT *seeds_idx, const IT &nb_seeds, FT *U, FT *DUr, FT *DUc, IT *V, IT label1 = 1)
    {
      IT idx;
      for (IT r = 0; r < nb_seeds; r++, label1++)
      {
	idx = seeds_idx[r];
	DUr[idx] = DUc[idx] = U[idx] = 0;
	_S[idx] = S_TRIAL;
	V[idx] = label1;
	_hp->push(idx);
      }
    }
    // ------------------------------------------------------------
    /**
     * @brief Add labelized grid nodes as boundary condition of the eikonal equation.
     * @param[in] seeds_idx Array of seeds represented by their indicies.
     * @param[in] nb_seeds Number of seeds in the array.
     * @param[in] U Distance (Minimal action) map.
     */
    template <class DelGraph>
    void addLabeledNodes(IT *seeds_idx, const IT &nb_seeds, FT *U, IT *V, DelGraph &G, IT label1 = 1)
    {
      IT idx;
      for (IT r = 0; r < nb_seeds; r++, label1++)
      {
	idx = seeds_idx[r];
	U[idx] = 0;
	_S[idx] = S_TRIAL;
	V[idx] = label1;
	_hp->push(idx);
	G.newNode(idx);
      }
    }
    // ------------------------------------------------------------
    inline void markNode(const IT &idx) { _S[idx] = S_MARK;  }
    // ------------------------------------------------------------
    /**
     * @brief Neighbor of a node in a simplex (i,j,n).
     * @param[in] index Index of the node i.
     * @param[in] radius Adjacency (0->4-adjacency, 1->8-adjacency).
     * @param[in] k Position of the node in the neighborhood of node j.
     * @param[in] l Flag (0 or 1) indicating one of the two possible other nodes of the simplices (i,j,n1) and (i,j,n2).
     */
    inline IT simplexNeighbor(const IT &index, const unsigned short &radius,
			      const unsigned short &k, unsigned short l)
    { return (index + _sn[radius][l][k]); }
    // ------------------------------------------------------------
    /**
     * @brief 4-Neighbor of a node in a simplex (i,j,n).
     * @param[in] index Index of the node i.
     * @param[in] k Position of the node in the neighborhood of node j.
     * @param[in] l Flag (0 or 1) indicating one of the two possible other nodes of the simplices (i,j,n1) and (i,j,n2).
     */
    inline IT simplexNeighbor4(const IT &index, const unsigned short &k, unsigned short l)
    { return (index + _sn4[l][k]); }
    // ------------------------------------------------------------
    /**
     * @brief 8-Neighbor of a node in a simplex (i,j,n).
     * @param[in] index Index of the node i.
     * @param[in] k Position of the node in the neighborhood of node j.
     * @param[in] l Flag (0 or 1) indicating one of the two possible other nodes of the simplices (i,j,n1) and (i,j,n2).
     */
    inline IT simplexNeighbor8(const IT &index, const unsigned short &k, unsigned short l)
    { return (index + _sn8[l][k]); }
    // ------------------------------------------------------------
  public:
    // ------------------------------------------------------------
    static const FT INFPLUS;
    // ------------------------------------------------------------
    /**
     * @brief Constructor.
     * @param[in] G Grid (must inherit from Grid8 or Grid4).
     */
    FM2(const GridClass &G)
      : _G(G), _S(NULL), S_TRIAL(-2), S_DEAD(2), _hp(NULL), _kvar(0), _P1(NULL), _P2(NULL), _medialgraph(NULL)
    {
      _sn4[0] = _sn4[1] = _sn8[0] = _sn8[1] = _wsn8[0] = _wsn8[1] = NULL;
      constructSimplexNeighbors();
      //constructSolutions();
    }
    // ------------------------------------------------------------
    /**
     * @brief Destructor.
     */
    ~FM2()
    {
      if (_hp) delete _hp;
      //if (_S) delete[] _S;
      if (_sn4[0]) delete[] _sn4[0]; if (_sn4[1]) delete[] _sn4[1];
      if (_sn8[0]) delete[] _sn8[0]; if (_sn8[1]) delete[] _sn8[1];
      if (_wsn8[0]) delete[] _wsn8[0]; if (_wsn8[1]) delete[] _wsn8[1];
      //if (_fct[0]) delete[] _fct[0]; if (_fct[1]) delete[] _fct[1];
      //if (_W[0]) delete[] _W[0]; if (_W[1]) delete[] _W[1];
      //if (_P1) delete[] _P1; if (_P2) delete[] _P2;
      if (_medialgraph) delete _medialgraph;
    }
    // ------------------------------------------------------------
    /**
     * @brief Geodesic distance from one source node s, and its gradient, such that U[s]=0.
     * @param[in] idx Index of the source node.
     * @param[in] pot Potential (scalar field) on the set of nodes. Negative values indicate nodes outside of the definition domain of U.
     * @param[out] U Geodesic distance from the source node.
     * @param[out] DUr Row derivative of the geodesic distance.
     * @param[out] DUc Column derivative of the geodesic distance.
     * @remark Boundary nodes are not treated by the FM (use method fillBoundary of Grid to modify the values if needed).
     * @par Example
     * @code
Grid4<long> g(nrows,ncols);
FM2<double,Grid4<long> > FM(g);
FM.distance(idx,P,U); // P is a potential
     * @endcode
     * Several distances U according to different potentials or domains.
     <table>
     <tr><td></td><td>@image html dgeoto1seed-pot.png "Potential P and source node"</td><td>@image html shapeto1seed.png "2D shape (P=-1 outside and P=1 inside)"</td></tr>
     <tr><td>@image html deucto1seed.png "Euclidean distance U (P=1)"</td><td>@image html dgeoto1seed.png "Geodesic distance U"</td><td>@image html dshapeto1seed.png "Shape inner distance U"</td></tr>
     </table>
     */
    void distance(IT idx, FT *pot, FT *U, FT *DUr = NULL, FT *DUc = NULL, IT *P1 = NULL, IT *P2 = NULL)
    {
      if (P1 && P2)
      {
	_P1 = P1; _P2 = P2;
	if (DUr && DUc)
	{
	  initPred(pot,U,DUr,DUc,_S);
	  addNodePred(idx,U,DUr,DUc);
	  propagatePred(pot,U,DUr,DUc);
	}
      }
      else
      {
	if (DUr && DUc)
	{
	  init(pot,U,DUr,DUc,_S);
	  addNode(idx,U,DUr,DUc);
	  propagate(pot,U,DUr,DUc);
	}
	else
	{
	  init(pot,U,_S);
	  addNode(idx,U);
	  propagate(pot,U);
	}
      }
      delete[] _S; _S = NULL;
    }
    // ------------------------------------------------------------
    void predecessorMap(IT idx, IT *P1, IT *P2, IT *P)
    {
      std::queue<IT> q;
      _G.fill(P,0);
      IT i1 = P1[idx], i2 = P2[idx];
      if (i1 > 0 && i1 != idx) { q.push(i1); P[i1] = 1; }
      if (i2 > 0 && i2 != idx) { q.push(i2); P[i2] = 1; }
      IT m = _G.size();
      while (!q.empty() && --m > 0)
      {
	idx = q.front();
	q.pop();
	if (P[idx] == 0) std::cerr << "err\n";
	i1 = P1[idx];
	i2 = P2[idx];
	if (i1 > 0 && i1 != idx && P[i1] == 0) { q.push(i1); P[i1] = 1; }
	if (i2 > 0 && i2 != idx && P[i2] == 0) { q.push(i2); P[i2] = 1; }
      }
    }
    // ------------------------------------------------------------
    /**
     * @brief Geodesic distance to a node set S, and its gradient, such that U[s]=0 for all s in S. Only one propagation is used to find the solution, which implies approximated values near the interface between fronts (Voronoi edges). Exact values are obtained using method distance.
     * @param[in] point_set_idx Array of source nodes (given by their indicies).
     * @param[in] nb_points Number of source nodes.
     * @param[in] pot Potential (scalar field). Negative values are considered outside of the definition domain of U.
     * @param[out] U Geodesic distance.
     * @param[out] DUr Row derivative of the geodesic distance.
     * @param[out] DUc Column derivative of the geodesic distance.
     */
    void distanceFast(IT *point_set_idx, IT nb_points, FT *pot, FT *U, FT *DUr = NULL, FT *DUc = NULL)
    {
      if (DUr)
      {
	init(pot,U,DUr,DUc,_S);
	addNodes(point_set_idx,nb_points,U,DUr,DUc);
	propagate(pot,U,DUr,DUc);
      }
      else
      {
	init(pot,U,_S);
	addNodes(point_set_idx,nb_points,U);
	propagate(pot,U);
      }
      delete[] _S; _S = NULL;
    }
    // ------------------------------------------------------------
    /**
     * @brief Geodesic distance to a node set S, and its gradient, such that U[s]=0 for all s in S. Several restricted propagations are used to compute the distance, one for each source node. A fast version can be obtained by using distanceFast, implying some approximations.
     * @param[in] point_set_idx Array of source nodes (given by their indicies).
     * @param[in] nb_points Number of source nodes.
     * @param[in] pot Potential (scalar field). Negative values are considered outside of the definition domain of U.
     * @param[out] U Geodesic distance.
     * @param[out] DUr Row derivative of the geodesic distance.
     * @param[out] DUc Column derivative of the geodesic distance.
     */
    void distance(IT *point_set_idx, IT nb_points, FT *pot, FT *U, FT *DUr = NULL, FT *DUc = NULL)
    {
      if (DUr)
      {
	// distance to first point
	init(pot,U,DUr,DUc,_S);
	addNode(point_set_idx[0],U,DUr,DUc);
	propagate(pot,U,DUr,DUc);
	// distance to other points computed iteratively
	FT *Upred = new FT[_G.size()];
	FT *DUrp =  new FT[_G.size()];
	FT *DUcp =  new FT[_G.size()];
	_hp->setData(Upred);
	for (IT i = 1; i < nb_points; i++)
	{
	  init(pot,Upred,DUrp,DUcp);
	  addNode(point_set_idx[i],Upred,DUrp,DUcp);
	  propagate(pot,U,DUr,DUc,Upred,DUrp,DUcp);
	}
	delete[] Upred; delete[] DUrp; delete[] DUcp;
      }
      else
      {
	// distance to first point
	init(pot,U,_S);
	addNode(point_set_idx[0],U);
	propagate(pot,U);
	// distance to other points computed iteratively
	FT *Upred = new FT[_G.size()];
	_hp->setData(Upred);
	for (IT i = 1; i < nb_points; i++)
	{
	  init(pot,Upred);
	  addNode(point_set_idx[i],Upred);
	  propagate(pot,U,Upred);
	}
	delete[] Upred;
      }
      delete[] _S; _S = NULL;
    }
    // ------------------------------------------------------------
    /**
     * @brief Update geodesic distance U by adding a new source node s with U[s]=0. To be used in iterative distance computation. The whole domain is not always explored during the propagation, which depends on the metric (potential), and on the values of the initial distance map U.
     * @param[in] idx Index of the source node.
     * @param[in] pot Potential (scalar field) on the set of nodes. Negative values indicate nodes outside of the definition domain of U.
     * @param[in,out] Uprev Geodesic distance already computed from some source nodes, with a distance method only (distanceFast should not be used).
     * @param[in,out] DUr Row derivative of Uprev.
     * @param[in,out] DUc Column derivative of Uprev.
     */
    void distanceUpdate(IT idx, FT *pot, FT *Uprev, FT *DUrp, FT *DUcp)
    {
      FT *U = new FT[_G.size()];
      FT *DUr = new FT[_G.size()];
      FT *DUc = new FT[_G.size()];
      init(pot,U,DUr,DUc,_S);
      addNode(idx,U,DUr,DUc);
      propagate(pot,Uprev,DUrp,DUcp,U,DUr,DUc);
      delete[] _S; _S = NULL;
      delete[] DUr; delete[] DUc; delete[] U;
    }
    // ------------------------------------------------------------
    /**
     * @brief Compute the Voronoi map, as well as the distance map and its gradient, in one propagation. The resulting Voronoi map contains errors on labels (as well as distances and gradients) near the interface between source nodes, that is Voronoi edges may not be well localized.
     * @param[in] point_set_idx Nodes represented by their indicies.
     * @param[in] nb_points Number of nodes.
     * @param[in] pot Potential map.
     * @param[out] V Voronoi map.
     * @param[out] U Distance (Minimal action) map.
     * @param[out] DUr Row gradient of the distance.
     * @param[out] DUc Column gradient of the distance.
     */
    void voronoiFast(IT *point_set_idx, IT nb_points, FT *pot, IT *V, FT *U, FT *DUr, FT *DUc)
    {
      _S = V;
      initLabelized(pot,U,DUr,DUc);
      addLabeledNodes(point_set_idx,nb_points,U,DUr,DUc);
      propagateLabeled(pot,U,DUr,DUc);
      _S = NULL;
    }
    // ------------------------------------------------------------
    /**
     * @brief Compute the Voronoi map, as well as the distance map and its gradient, in one propagation. The resulting Voronoi map contains errors on labels (as well as distances and gradients) near the interface between source nodes, that is Voronoi edges may not be well localized.
     * @param[in] point_set_idx Nodes represented by their indicies.
     * @param[in] nb_points Number of nodes.
     * @param[in] pot Potential map.
     * @param[out] V Voronoi map.
     * @param[out] U Distance (Minimal action) map.
     * @param[out] DUr Row gradient of the distance.
     * @param[out] DUc Column gradient of the distance.
     */
    void voronoiFast_old(IT *point_set_idx, IT nb_points, FT *pot, IT *V, FT *U, FT *DUr, FT *DUc)
    {
      init(pot,U,DUr,DUc,_S);
      addLabeledNodes(point_set_idx,nb_points,U,DUr,DUc,V);
      propagate(pot,U,DUr,DUc,V);
      delete[] _S; _S = NULL;
    }
    // ------------------------------------------------------------
    /**
     * @brief Geodesic distance to a node set S, and its gradient, such that U[s]=0 for all s in S. Several restricted propagations are used to compute the distance, one for each source node. A fast version can be obtained by using distanceFast, implying some approximations. Provides the same result as voronoiSlow but with restricted domain exploration during propagations.
     * @param[in] point_set_idx Array of source nodes (given by their indicies).
     * @param[in] nb_points Number of source nodes.
     * @param[in] pot Potential (scalar field). Negative values are considered outside of the definition domain of U.
     * @param[out] U Geodesic distance.
     * @param[out] DUr Row derivative of the geodesic distance.
     * @param[out] DUc Column derivative of the geodesic distance.
     */
    void voronoi(IT *point_set_idx, IT nb_points, FT *pot, IT *V, FT *U, FT *DUr, FT *DUc)
    {
      // distance to first point
      IT lab_idx = 1;
      init(pot,U,DUr,DUc,_S,V,lab_idx);
      addNode(point_set_idx[0],U,DUr,DUc);
      propagate(pot,U,DUr,DUc);
      // distance to other points computed iteratively
      FT *Upred = new FT[_G.size()];
      FT *DUrp =  new FT[_G.size()];
      FT *DUcp =  new FT[_G.size()];
      _hp->setData(Upred);
      for (IT i = 1; i < nb_points; i++)
      {
	lab_idx++;
	init(pot,Upred,DUrp,DUcp);
	addNode(point_set_idx[i],Upred,DUrp,DUcp);
	propagate(pot,V,lab_idx,U,DUr,DUc,Upred,DUrp,DUcp);
      }
      //_G.copy(_S,V);
      delete[] Upred; delete[] DUrp; delete[] DUcp;
      delete[] _S; _S = NULL;
    }
    // ------------------------------------------------------------
    /**
     * @brief Geodesic distance to a node set S, and its gradient, such that U[s]=0 for all s in S. Several propagations are used to compute the distance, one for each source node.
     * @param[in] point_set_idx Array of source nodes (given by their indicies).
     * @param[in] nb_points Number of source nodes.
     * @param[in] pot Potential (scalar field). Negative values are considered outside of the definition domain of U.
     * @param[out] U Geodesic distance.
     * @param[out] DUr Row derivative of the geodesic distance.
     * @param[out] DUc Column derivative of the geodesic distance.
     */
    void voronoiSlow(IT *point_set_idx, IT nb_points, FT *pot, IT *V, FT *U, FT *DUr, FT *DUc)
    {
      // distance to first point
      IT lab_idx = 1, idx;
      init(pot,U,DUr,DUc,_S,V,lab_idx);
      addNode(point_set_idx[0],U,DUr,DUc);
      propagate(pot,U,DUr,DUc);
      // distance to other points computed iteratively
      FT *Upred = new FT[_G.size()];
      FT *DUrp =  new FT[_G.size()];
      FT *DUcp =  new FT[_G.size()];
      _hp->setData(Upred);
      for (IT i = 1; i < nb_points; i++)
      {
	lab_idx++;
	init(pot,Upred,DUrp,DUcp);
	addNode(point_set_idx[i],Upred,DUrp,DUcp);
	propagate(pot,Upred,DUrp,DUcp);
	for_grid_interior(_G,idx,IT)
	{
	  if (Upred[idx] < U[idx]) { U[idx] = Upred[idx]; V[idx] = lab_idx; }
	}
      }
      delete[] Upred; delete[] DUrp; delete[] DUcp;
      delete[] _S; _S = NULL;
    }
    // ------------------------------------------------------------
    /**
     * @brief Update geodesic Voronoi map V, with distance and its gradient, by adding a new source node s with U[s]=0. To be used in iterative distance computation. The whole domain is not always explored during the propagation, which depends on the metric (potential), and on the values of the initial distance map U.
     * @param[in] idx Index of the source node.
     * @param[in] lab_idx Label of the source node.
     * @param[in] pot Potential (scalar field) on the set of nodes. Negative values indicate nodes outside of the definition domain of U.
     * @param[in] Vprev Voronoi map already computed.
     * @param[in,out] Uprev Geodesic distance already computed from some source nodes, with a distance method only (distanceFast should not be used).
     * @param[in,out] DUr Row derivative of Uprev.
     * @param[in,out] DUc Column derivative of Uprev.
     */
    void voronoiUpdate(IT idx, IT lab_idx, FT *pot, IT *Vprev, FT *Uprev, FT *DUrp, FT *DUcp)
    {
      FT *U = new FT[_G.size()];
      FT *DUr = new FT[_G.size()];
      FT *DUc = new FT[_G.size()];
      init(pot,U,DUr,DUc,_S);
      addNode(idx,U,DUr,DUc);
      propagate(pot,Vprev,lab_idx,Uprev,DUrp,DUcp,U,DUr,DUc);
      delete[] _S; _S = NULL;
      delete[] DUr; delete[] DUc; delete[] U;
    }
    // ------------------------------------------------------------
    /*
     * @brief Compute the Voronoi map, as well as the distance map and its gradient, in one propagation. The resulting Voronoi map contains errors on labels (as well as distances and gradients) near the interface between nodes, that is Voronoi edges may not be well localized.
     * @param[in] point_set_idx Nodes represented by their indicies.
     * @param[in] nb_points Number of nodes.
     * @param[in] pot Potential map.
     * @param[out] V Voronoi map.
     * @param[out] U Distance (Minimal action) map.
     * @param[out] DUr Row gradient of the distance.
     * @param[out] DUc Column gradient of the distance.

    template <class DelGraph>
    void voronoiFast(IT *point_set_idx, IT nb_points, FT *pot, IT *V,  DelGraph &DG, FT *U, FT *DUr, FT *DUc)
    {
      constructDots();
      init(pot,U,DUr,DUc,_S);
      DG.newNode(0);
      addLabeledNodes(point_set_idx,nb_points,U,V,DG);
      propagate(pot,U,DUr,DUc,V,DG);
      delete[] _S; _S = NULL;
    }*/
    // ------------------------------------------------------------
    /**
     * @brief Compute the Voronoi map, as well as the distance map and its gradient, in one propagation. The resulting Voronoi map contains errors on labels (as well as distances and gradients) near the interface between nodes, that is Voronoi edges may not be well localized.
     * @param[in] point_set_idx Nodes represented by their indicies.
     * @param[in] nb_points Number of nodes.
     * @param[in] pot Potential map.
     * @param[out] V Voronoi map.
     * @param[out] DG Delaunay graph.
     * @param[out] U Distance (Minimal action) map.
     * @param[out] DUr Row gradient of the distance.
     * @param[out] DUc Column gradient of the distance.
     */
    template <class DelGraph>
    void voronoiFast(IT *point_set_idx, IT nb_points, FT *pot, IT *V,  DelGraph &DG, FT *U, FT *DUr, FT *DUc)
    {
      _S = V;
      DG.reserve(nb_points+2); // labels begin at 2 (0 is the border, 1 is unused)
      initLabelized(pot,U,DUr,DUc);
      addLabeledNodes(point_set_idx,nb_points,U,DUr,DUc);
      propagateLabeled(pot,U,DUr,DUc,DG);
      // std::vector<typename DelGraph::Edge*> le;
      // DG.edges(le);
      // _G.fillInterior(V,0);
      // for (typename std::vector<typename DelGraph::Edge*>::iterator eit = le.begin(), end = le.end(); eit != end; ++eit)
      // {
      // 	// IT r,c;
      // 	// _G.point((*eit)->idx1,r,c);
      //  	// std::cerr << (*eit)->lab1 << " " << (*eit)->lab2 << " " << (*eit)->d << " " << r+1 << " " << c+1 << "\n";
      // 	// _G.point((*eit)->idx2,r,c);
      // 	// std::cerr << r+1 << " " << c+1 << "\n";
      // 	(*eit)->sort();
      // 	typename DelGraph::Edge* e = *eit;
      // 	IT lab = 1;
      // 	const typename DelGraph::Edge::VPairs &ep = e->vpairs;
      // 	for (IT i = 0; i < ep.size(); i++)
      // 	{
      // 	  if (V[ep[i].idx1] == 0)
      // 	  {
      // 	    if (V[ep[i].idx2] == 0)
      // 	    {
      // 	      bool t = true;
      // 	      IT labtmp = -1;
      // 	      for (unsigned short k = 0; k < 8; k++)
      // 	      {
      // 		IT idxn = _G.neighbor8(ep[i].idx1,k);
      // 		if (V[idxn] != 0)
      // 		{
      // 		  if (labtmp < 0) labtmp = V[idxn];
      // 		  else
      // 		    if (labtmp != V[idxn]) std::cerr << labtmp << " " << V[idxn] << " detect\n";
      // 		  V[ep[i].idx1] = V[idxn]; V[ep[i].idx2] = V[idxn]; t = false;
      // 		}
      // 	      }
      // 	      for (unsigned short k = 0; k < 8; k++)
      // 	      {
      // 		IT idxn = _G.neighbor8(ep[i].idx2,k);
      // 		if (V[idxn] != 0) //{ V[(*eit)->vpairs[i].idx1] = V[idxn]; V[(*eit)->vpairs[i].idx2] = V[idxn]; t = false; }
      // 		{
      // 		  if (labtmp < 0) labtmp = V[idxn];
      // 		  else
      // 		    if (labtmp != V[idxn]) std::cerr << labtmp << " " << V[idxn] << " detect\n";
      // 		  V[ep[i].idx1] = V[idxn]; V[ep[i].idx2] = V[idxn]; t = false;
      // 		}
      // 	      }
      // 	      if (t)
      // 	      {
      // 		V[ep[i].idx1] = lab;
      // 		V[ep[i].idx2] = lab;
      // 		lab++;
      // 	      }
      // 	    }
      // 	    else V[ep[i].idx1] = V[ep[i].idx2];
      // 	  }
      // 	  else
      // 	  {
      // 	    if (V[ep[i].idx2] == 0) V[ep[i].idx2] = V[ep[i].idx1];
      // 	    else
      // 	    {

      // 	    }
      // 	  }
      // 	  // ------------------------------------------------------------
      // 	  if (V[ep[ep.size()-i-1].idx1] == 0)
      // 	  {
      // 	    if (V[ep[ep.size()-i-1].idx2] == 0)
      // 	    {
      // 	      bool t = true;
      // 	      IT labtmp = -1;
      // 	      for (unsigned short k = 0; k < 8; k++)
      // 	      {
      // 		IT idxn = _G.neighbor8(ep[ep.size()-i-1].idx1,k);
      // 		if (V[idxn] != 0)
      // 		{
      // 		  if (labtmp < 0) labtmp = V[idxn];
      // 		  else
      // 		    if (labtmp != V[idxn]) std::cerr << labtmp << " " << V[idxn] << " detect\n";
      // 		  V[ep[ep.size()-i-1].idx1] = V[idxn]; V[ep[ep.size()-i-1].idx2] = V[idxn]; t = false;
      // 		}
      // 	      }
      // 	      for (unsigned short k = 0; k < 8; k++)
      // 	      {
      // 		IT idxn = _G.neighbor8(ep[ep.size()-i-1].idx2,k);
      // 		if (V[idxn] != 0) //{ V[(*eit)->vpairs[i].idx1] = V[idxn]; V[(*eit)->vpairs[i].idx2] = V[idxn]; t = false; }
      // 		{
      // 		  if (labtmp < 0) labtmp = V[idxn];
      // 		  else
      // 		    if (labtmp != V[idxn]) std::cerr << labtmp << " " << V[idxn] << " detect\n";
      // 		  V[ep[ep.size()-i-1].idx1] = V[idxn]; V[ep[ep.size()-i-1].idx2] = V[idxn]; t = false;
      // 		}
      // 	      }
      // 	      if (t)
      // 	      {
      // 		V[ep[ep.size()-i-1].idx1] = lab;
      // 		V[ep[ep.size()-i-1].idx2] = lab;
      // 		lab++;
      // 	      }
      // 	    }
      // 	    else V[ep[ep.size()-i-1].idx1] = V[ep[ep.size()-i-1].idx2];
      // 	  }
      // 	  else
      // 	  {
      // 	    if (V[ep[ep.size()-i-1].idx2] == 0) V[ep[ep.size()-i-1].idx2] = V[ep[ep.size()-i-1].idx1];
      // 	    else
      // 	    {

      // 	    }
      // 	  }
      // 	}
      // }
    }
    // ------------------------------------------------------------
    /**
     * @brief Compute the Voronoi map, as well as the distance map and its gradient, in several restricted propagations (one for each source node).
     * @param[in] point_set_idx Nodes represented by their indicies.
     * @param[in] nb_points Number of nodes.
     * @param[in] pot Potential map.
     * @param[out] V Voronoi map.
     * @param[out] DG Delaunay graph.
     * @param[out] U Distance (Minimal action) map.
     * @param[out] DUr Row gradient of the distance.
     * @param[out] DUc Column gradient of the distance.
     */
    template <class DelGraph>
    void voronoi(IT *point_set_idx, IT nb_points, FT *pot, IT *V,  DelGraph &DG, FT *U, FT *DUr, FT *DUc)
    {
      // distance to first point
      DG.reserve(nb_points+2); // labels begin at 2 (0 is the border, 1 is unused)
      IT lab_idx = 2;
      init(pot,U,DUr,DUc,_S,V,lab_idx);
      addNode(point_set_idx[0],U,DUr,DUc);
      propagate(pot,U,DUr,DUc);
      // distance to other points computed iteratively
      FT *Upred = new FT[_G.size()];
      FT *DUrp =  new FT[_G.size()];
      FT *DUcp =  new FT[_G.size()];
      _hp->setData(Upred);
      for (IT i = 1; i < nb_points; i++)
      {
	lab_idx++;
	init(pot,Upred,DUrp,DUcp);
	addNode(point_set_idx[i],Upred,DUrp,DUcp);
	propagate(pot,V,lab_idx,U,DUr,DUc,Upred,DUrp,DUcp,DG);
      }
      delete[] Upred; delete[] DUrp; delete[] DUcp;
      delete[] _S; _S = NULL;
      std::vector<typename DelGraph::Edge*> le;
      DG.edges(le);
      for (typename std::vector<typename DelGraph::Edge*>::iterator eit = le.begin(), end = le.end(); eit != end; ++eit)
      {
      	IT r,c;
      	_G.point((*eit)->idx1,r,c);
       	std::cerr << (*eit)->lab1 << " " << (*eit)->lab2 << " " << (*eit)->d << " " << r+1 << " " << c+1 << "\n";
      	_G.point((*eit)->idx2,r,c);
      	std::cerr << r+1 << " " << c+1 << "\n";
      }
    }
    // ------------------------------------------------------------
    /**
     * @brief Compute the Voronoi map, as well as the distance map and its gradient, in one propagation. The resulting Voronoi map contains errors on labels (as well as distances and gradients) near the interface between nodes, that is Voronoi edges may not be well localized.
     * @param[in] point_set_idx Nodes represented by their indicies.
     * @param[in] nb_points Number of nodes.
     * @param[in] pot Potential map.
     * @param[out] V Voronoi map.
     * @param[out] U Distance (Minimal action) map.
     * @param[out] DUr Row gradient of the distance.
     * @param[out] DUc Column gradient of the distance.
     */
    template <class DelGraph>
    void restrictedVoronoiFast(IT *point_set_idx, IT nb_points, FT *pot, DelGraph &DG, FT *U, FT *DUr, FT *DUc, IT *V)
    {
      _S = V;
      initLabelized(pot,U,DUr,DUc);
      addLabeledNodes(point_set_idx,nb_points,U);
      propagateLabeled(pot,U,DUr,DUc,DG,nb_points);

      initLabelized(pot,U,DUr,DUc);
      std::vector<typename DelGraph::Edge*> le;
      DG.edges(le);
      for (typename std::vector<typename DelGraph::Edge*>::iterator eit = le.begin(), end = le.end(); eit != end; ++eit)
      {
       	std::cerr << (*eit)->idx1 << " " << (*eit)->idx2 << "\n";
       	addLabeledNodes((*eit)->idx1,(*eit)->idx2,U);
      }
      DelGraph DG2(nb_points+2+le.size());
      propagateLabeled(pot,U,DUr,DUc,DG2,le.size());
    }
    // ------------------------------------------------------------
    /**
     * @brief Geodesic distance from one source node s to a target node, and its gradient, with the Fast Marching algorithm and such that U[s]=0.
     * @param[in] idx Index of the source node.
     * @param[in] idx_to Index of the target node.
     * @param[in] pot Potential (scalar field) on the set of nodes. Negative values indicate the outside of the definition domain of U.
     * @param[out] U Geodesic distance to the source node.
     * @param[out] DUr Row derivative of the geodesic distance to the source node.
     * @param[out] DUc Column derivative of the geodesic distance to the source node.
     * @return The distance from idx to idx_to.
     */
    FT distanceTo(IT idx, IT idx_to, FT *pot, FT *U, FT *DUr, FT *DUc)
    {
      init(pot,U,_S);
      addNode(idx,U);
      FT res = propagate(idx_to,pot,U,DUr,DUc);
      delete[] _S; _S = NULL;
      return res;
    }
    // ------------------------------------------------------------
    /**
     * @brief Geodesic distance to one source node s, and its gradient, with the Fast Marching algorithm and such that U[s]=0.
     * @param[in] idx Index of the source node.
     * @param[in] pot Potential (scalar field) on the set of nodes. Negative values indicate the outside of the definition domain of U.
     * @param[out] U Geodesic distance to the source node.
     * @param[out] DUr Row derivative of the geodesic distance to the source node.
     * @param[out] DUc Column derivative of the geodesic distance to the source node.
     */
    template <class WeightClass>
    void directedDistance(IT idx, const WeightClass &W, FT *U, FT *DUr, FT *DUc)
    {
      constructWSimplex8Neighbors();
      init(U,_S);
      addNode(idx,U,DUr,DUc);
      directedPropagation(W,U,DUr,DUc);
      delete[] _S; _S = NULL;
    }
    // ------------------------------------------------------------
    void medialSet(IT idx, IT lab_idx, FT *pot, FT *Uprev, FT *DUr, FT *DUc, IT *V, FT *Unew, FT *DUrnew, FT *DUcnew, FT *DUS)
    {
      _G.copy(Uprev,Unew);
      _G.copy(DUr,DUrnew);
      _G.copy(DUc,DUcnew);
      _G.fill(DUS,(FT)0);
      constructDots();
      _label = lab_idx;
      _medialgraph = new IndexedCMap<IT,Node<IT,FT> >;
      FT *U = new FT[_G.size()];
      FT *DUrt = new FT[_G.size()];
      FT *DUct = new FT[_G.size()];
      // _P1 = new IT[_G.size()]; _P2 = new IT[_G.size()];
      init(pot,U,_S);
      addNode(idx,U);
      _nbAcceptedInside = 0;
      while (!_hp->empty())
      {
	idx = _hp->pop();
	if (Unew[idx] < U[idx]) // not inside
	{
	  _S[idx] = S_DEAD_OUTSIDE;
	  // updateNeighborsOutside<GridClass>(idx,pot,U,DUrt,DUct,Uprev,DUr,DUc,V,DUS);
	}
	else  // inside
	{
	  _S[idx] = S_DEAD_INSIDE;
	  V[idx] = _label;
	  _nbAcceptedInside++;
	  // updateNeighborsInside<GridClass>(idx,pot,U,DUrt,DUct,V);
	  Unew[idx] = U[idx]; DUrnew[idx] = DUrt[idx]; DUcnew[idx] = DUct[idx];
	}
	if (_nbAcceptedInside <= 0) {break;}  // all accepted inside have been computed
      }

      // medial graph
      std::vector<std::vector<NodeDt<IT,FT> > > v;
      _medialgraph->connectedComponents(v,_G);
      // for (IT j = 0; j < v.size(); j++)
      // {
      // 	IT r,c;
      // 	for (IT i = 0; i < v[j].size(); i++)
      // 	{
      // 	  _G.point(v[j][i],r,c);
      // 	  std::cerr << "[" << r+1.5 << "," << c+1.5 << "];";
      // 	}
      // 	std::cerr << "\n";
      // }
      // std::vector<IT> c;
      // _medialgraph->degree1Nodes(c);
      // std::cerr << c.size() << " nodes of degree 1\n";
      // // for each degree 1 node, try to find a 8-neighbor to close the contour
      // for (IT i = 0; i < c.size(); i++)
      // {
      // 	FOREACH_8DIAG_NEIGHBOR(_kvar)
      // 	{
      // 	  _idx_n = _G.neighbor8(c[i],_kvar);
      // 	  if (DUS[_idx_n] >= 0)
      // 	    if (_medialgraph->node(c[i])->degree == 1) _medialgraph->link(c[i],_idx_n);
      // 	}
      // }
      // c.clear();
      // _medialgraph->degree1Nodes(c);
      // std::cerr << c.size() << " nodes of degree 1\n";

      // Node<IT,FT>* nd = _medialgraph->node();
      // if (nd)
      // {
      // 	IT didx = nd->dart;
      // 	if (didx != -1)
      // 	{
      // 	  Dart<IT>* d = _medialgraph->dart(didx);
      // 	  IT a = d->alpha;
      // 	  if (a != -1)
      // 	  {
      // 	    d = _medialgraph->dart(a);
      // 	    std::cerr << d->node << "\n";
      // 	  }
      // 	}
      // }

      delete[] U; delete[] _S; _S = NULL; delete[] DUrt; delete[] DUct;
      //delete[] _P1; delete[] _P2; _P1 = NULL; _P2 = NULL;
      delete _medialgraph; _medialgraph = NULL;
    }
    // ------------------------------------------------------------
    /*
     * @brief Compute the geodesic distance to one source node s (approximated solution of the Eikonal equation) with the Fast Marching algorithm and such that U[s]=0.
     * @param[in] idx Index of the source node.
     * @param[in] pot Potential (scalar) field on the set of nodes. Negative values indicate the outside of the definition domain of U.
     * @param[out] U Geodesic distance to the source node (arrival time).
     */
    /*void constrainedDistance2NN(IT idx, IT idxn1, IT idxn2, FT *pot, FT *U, FT *DUr, FT *DUc)
    {
      init(pot,U,_S); addNode(idx,U);
      IT k = 0;
      while (!_hp->empty())
      {
	idx = _hp->pop();
	if (idx == idxn1 || idx == idxn2) { k++; if (k == 2) break; }
	_S[idx] = S_DEAD;
	updateNeighbors<GridClass>(idx,pot,U,DUr,DUc);
      }
      while (!_hp->empty()) { idx = _hp->pop(); _U[idx] = INFPLUS; }
      delete[] _S; _S = NULL;
      }*/
    // ------------------------------------------------------------
    /**
     * @brief Farthest point sampling of the grid domain, according to the potential.
     * @param[in] pot Potential (scalar) field on the set of nodes. Negative values indicate the outside of the definition domain of U.
     * @param[out] U Geodesic distance to the source node (arrival time).
     */
    template <class NodeIdxContainer>
    void fps(FT *pot, IT nb_points, NodeIdxContainer &point_set_idx, FT *U, FT *DUr, FT *DUc, IT *V)
    {
      // distance to first point
      IT lab_idx = 1, idx = _G.firstInteriorRowFirstNode();
      init(pot,U,DUr,DUc,_S,V,lab_idx);
      addNode(idx,U,DUr,DUc);
      point_set_idx.push_back(idx);
      propagate(pot,U,DUr,DUc);
      // distance to the 3 other corner points, computed iteratively
      FT *Upred = new FT[_G.size()];
      FT *DUrp =  new FT[_G.size()];
      FT *DUcp =  new FT[_G.size()];
      _hp->setData(Upred);
      point_set_idx.push_back(_G.firstInteriorRowLastNode());
      point_set_idx.push_back(_G.lastInteriorRowFirstNode());
      point_set_idx.push_back(_G.lastInteriorRowLastNode());
      typename NodeIdxContainer::iterator nit = point_set_idx.begin(), end = point_set_idx.end();
      for (; nit != end; ++nit)
      {
	lab_idx++;
	init(pot,Upred,DUrp,DUcp);
	addNode(*nit,Upred,DUrp,DUcp);
	propagate(pot,V,lab_idx,U,DUr,DUc,Upred,DUrp,DUcp);
      }
      // distance to farthest points, computed iteratively
      nb_points -= 4;
      while (nb_points > 0)
      {
	lab_idx++;
	initFPS(pot,Upred,DUrp,DUcp,U,idx); // also finds the farthest point of previously computed U
	addNode(idx,Upred,DUrp,DUcp);
	point_set_idx.push_back(idx);
	propagate(pot,V,lab_idx,U,DUr,DUc,Upred,DUrp,DUcp);
	nb_points--;
      }
      delete[] Upred; delete[] DUrp; delete[] DUcp;
      delete[] _S; _S = NULL;
    }
    // ------------------------------------------------------------
    /**
     * @brief Path starting at idx_from and ending at the nearest node satisfying Pred[n]=n (U[n]=0).
     * @param[in] idx_from Index of the starting node.
     * @param[in] DUr Row gradient of distance array.
     * @param[in] DUc Column gradient of distance array.
     * @param[out] path Minimal path (Container must have a method push_back(elt)).
     * @param[in] terminal_nodes Add terminal nodes to the path (by default) or not.
     * @par Example
     * @image html dusum.png "Norm of DU sum"
     */
    template <class Container>
    void minimalPath(IT idx_from, FT *DUr, FT *DUc, Container &path,
		     bool terminal_nodes = true, FT dt = 0.6)
    { _G.odeEuler(idx_from,DUr,DUc,path,terminal_nodes,dt); }
    // ------------------------------------------------------------
  };
  // ===================================================
  // static members declaration
  // ===================================================
  template <typename FT, class GrdClass>
  const FT FM2<FT,GrdClass>::_SQRT2 = std::sqrt((FT)2.);

  template <typename FT, class GrdClass>
  const FT FM2<FT,GrdClass>::SIGNR4[4][2] = {{-1.0,1.0},{1.0,1.0},{-1.0,1.0},{-1.0,-1.0}};

  template <typename FT, class GrdClass>
  const FT FM2<FT,GrdClass>::SIGNC4[4][2] = {{1.0,1.0},{-1.0,1.0},{-1.0,-1.0},{-1.0,1.0}};

  template <typename FT, class GrdClass>
  const FT FM2<FT,GrdClass>::SIGNR8[8][2] = {{-1.0,1.0},{1.0,1.0},{1.0,1.0},{1.0,1.0},{-1.0,1.0},{-1.0,-1.0},{-1.0,-1.0},{-1.0,-1.0}};

  template <typename FT, class GrdClass>
  const FT FM2<FT,GrdClass>::SIGNC8[8][2] = {{1.0,1.0},{1.0,1.0},{-1.0,1.0},{-1.0,-1.0},{-1.0,-1.0},{-1.0,-1.0},{-1.0,1.0},{1.0,1.0}};

  template <typename FT, class GrdClass>
  const typename GrdClass::IndexType FM2<FT,GrdClass>::S_FAR = -1;

  template <typename FT, class GrdClass>
  const typename GrdClass::IndexType FM2<FT,GrdClass>::S_BORDER = 0;

  template <typename FT, class GrdClass>
  const typename GrdClass::IndexType FM2<FT,GrdClass>::S_MARK = 1;

  template <typename FT, class GrdClass>
  const typename GrdClass::IndexType FM2<FT,GrdClass>::S_DEAD_INSIDE = 3;

  template <typename FT, class GrdClass>
  const typename GrdClass::IndexType FM2<FT,GrdClass>::S_DEAD_OUTSIDE = 4;

  template <typename FT, class GrdClass>
  const FT FM2<FT,GrdClass>::INFPLUS = std::numeric_limits<FT>::max();
}

#endif
