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
 * @file cmap.hh
 */
// =========================================================================
#ifndef __CMAP_HH__
#define __CMAP_HH__

#include <map>
#include <vector>
#include <limits>

namespace GeoComp
{
  // ------------------------------------------------------------
  // ------------------------------------------------------------
  template <typename IT = int>
  class HalfEdge
  {
  public:
    // ------------------------------------------------------------
    HalfEdge *opposite;
    HalfEdge *next;
    IT nodeIdx;
    // ------------------------------------------------------------
    HalfEdge() : opposite(NULL), next(this), nodeIdx(-1) { }
    // ------------------------------------------------------------
    HalfEdge(const IT &nidx) : opposite(NULL), next(this), nodeIdx(nidx) { }
    // ------------------------------------------------------------
    HalfEdge(const IT &nidx, HalfEdge *op) : opposite(op), next(this), nodeIdx(nidx) { }
    // ------------------------------------------------------------
    HalfEdge(const IT &nidx, HalfEdge *op, HalfEdge *n) : opposite(op), next(n), nodeIdx(nidx) { }
    // ------------------------------------------------------------
    ~HalfEdge() { }
  };
  
  // ------------------------------------------------------------
  template <typename IT, class HEDataCl>
  class HalfEdgeData
  {
  public:
    typedef HEDataCl HEDataClass;
    // ------------------------------------------------------------
    HalfEdgeData *opposite;
    HalfEdgeData *next;
    IT nodeIdx;
    // ------------------------------------------------------------
    HEDataClass *data;
    // ------------------------------------------------------------
    HalfEdgeData() : opposite(NULL), next(this), nodeIdx(-1), data(NULL) { }
    // ------------------------------------------------------------
    HalfEdgeData(const IT &nidx)
      : opposite(NULL), next(this), nodeIdx(nidx), data(NULL) { }
    // ------------------------------------------------------------
    HalfEdgeData(const IT &nidx, HalfEdgeData<IT,HEDataClass> *op)
      : opposite(op), next(this), nodeIdx(nidx), data(NULL) { }
    // ------------------------------------------------------------
    HalfEdgeData(const IT &nidx, HalfEdgeData<IT,HEDataClass> *op, HalfEdgeData<IT,HEDataClass> *n)
      : opposite(op), next(n), nodeIdx(nidx), data(NULL) { }
    // ------------------------------------------------------------
    ~HalfEdgeData()
    {
      if (data)
      {
	if (opposite && opposite->data == data) opposite->data = NULL;
	delete data;
      }
    }
  };

  // ------------------------------------------------------------
  // ------------------------------------------------------------
  template <typename IT, class HECl>
  class NodeBase
  {
  public:
    typedef IT IdxType;
    typedef HECl HEClass;
    HEClass *halfEdge;  // incident HE
    NodeBase() : halfEdge(NULL) { }
    ~NodeBase() { }
  };

  // ------------------------------------------------------------
  // ------------------------------------------------------------
  template <typename IT, class HECl, class NodeDataClass>
  class NodeBaseData
  {
  public:
    typedef NodeDataClass NodeData;
    // ------------------------------------------------------------
    typedef IT IdxType;
    typedef HECl HEClass;
    HEClass *halfEdge;  // incident HE
    NodeData *data;
    NodeBaseData() : halfEdge(NULL), data(NULL) { }
    ~NodeBaseData()
    {
      if (data) delete data;
    }
  };
  
  // ------------------------------------------------------------
  // ------------------------------------------------------------
  template <class NdCl, class NdContainer = std::vector<NdCl*> >
  class HEGraph
  {
  public:
    typedef NdCl NodeClass;
    typedef typename NodeClass::HEClass HEClass;
    typedef typename NodeClass::IdxType IT;
    typedef typename std::vector<HEClass*>::iterator HEIterator;
    typedef typename NdContainer::iterator NodeIterator;
    // ------------------------------------------------------------
  protected:
    NdContainer _nodes;
    std::vector<HEClass*> _hes;
  public:
    // ------------------------------------------------------------
    HEGraph() { }
    // ------------------------------------------------------------
    ~HEGraph()
    {
      for (NodeIterator it = _nodes.begin(), end = _nodes.end(); it != end; ++it) delete *it;
      for (HEIterator it = _hes.begin(), end = _hes.end(); it != end; ++it) delete *it;
    }
    // ------------------------------------------------------------
    inline HEIterator halfEdgesBegin() { return _hes.begin(); }
    // ------------------------------------------------------------
    inline HEIterator halfEdgesEnd() { return _hes.end(); }
    // ------------------------------------------------------------
    inline NodeClass* isNode(const IT &idx) { return (idx < _nodes.size() ? _nodes[idx] : NULL); }
    // ------------------------------------------------------------
    inline NodeClass* node(const IT &idx) { return _nodes[idx]; }
    // ------------------------------------------------------------
    inline HEClass* incidentHalfEdge(const IT &idx) { return _nodes[idx]->halfEdge; }
    // ------------------------------------------------------------
    inline void newNode(NodeClass* n) { _nodes.push_back(n); }
    // ------------------------------------------------------------
    NodeClass* newNode()
    {
      NodeClass *res = new NodeClass;
      newNode(res);
      return res;
    }
    // ------------------------------------------------------------
    HEClass* newEdge(const IT &n1_idx, const IT &n2_idx)
    {
      // create halfedges
      HEClass *d12 = new HEClass(n1_idx);
      HEClass *d21 = new HEClass(n2_idx,d12);
      d12->opposite = d21;
      // insert them into array of halfedges
      _hes.push_back(d12);
      _hes.push_back(d21);
      // update orbit of n1
      NodeClass *n = node(n1_idx);
      HEClass *dn = n->halfEdge;
      n->halfEdge = d12;
      if (dn)
      {
	d12->next = dn->next;
	dn->next = d12;
      }
      // update orbit of n2
      n = node(n2_idx);
      dn = n->halfEdge;
      n->halfEdge = d21;
      if (dn)
      {
	d21->next = dn->next;
	dn->next = d21;
      }
      return d12;
    }
    // ------------------------------------------------------------
    HEClass* isNeighbor(const IT &n1_idx, const IT &n2_idx)
    {
      HEClass *he = incidentHalfEdge(n1_idx);
      if (he == NULL) return NULL;
      HEClass *end = he;
      do
      {
	if (he->opposite->nodeIdx == n2_idx) return he;
	he = he->next;
      } while (he != end);
      return NULL;
    }
    // ------------------------------------------------------------
  };

  // ------------------------------------------------------------
  template <class NdCl, class NdContainer = std::map<typename NdCl::IdxType,NdCl*> >
  class IndexedHEGraph
  {
  public:
    typedef NdCl NodeClass;
    typedef typename NodeClass::HEClass HEClass;
    typedef typename NodeClass::IdxType IT;
    typedef typename std::vector<HEClass*>::iterator HEIterator;
    typedef typename NdContainer::iterator NodeIterator;
    // ------------------------------------------------------------
  protected:
    NdContainer _nodes;
    std::vector<HEClass*> _hes;
  public:
    // ------------------------------------------------------------
    IndexedHEGraph() { }
    // ------------------------------------------------------------
    ~IndexedHEGraph()
    {
      for (NodeIterator it = nodesBegin(), end = nodesEnd(); it != end; ++it) delete it->second;
      for (HEIterator it = _hes.begin(), end = _hes.end(); it != end; ++it) delete *it;
    }
    // ------------------------------------------------------------
    inline HEIterator halfEdgesBegin() { return _hes.begin(); }
    // ------------------------------------------------------------
    inline HEIterator halfEdgesEnd() { return _hes.end(); }
    // ------------------------------------------------------------
    inline NodeIterator nodesBegin() { return _nodes.begin(); }
    // ------------------------------------------------------------
    inline NodeIterator nodesEnd() { return _nodes.end(); }
    // ------------------------------------------------------------
    inline NodeClass* isNode(const IT &idx)
    {
      NodeIterator nit = _nodes.find(idx);
      return (nit == _nodes.end() ? NULL : nit->second);
    }
    // ------------------------------------------------------------
    inline NodeClass* node(const IT &idx) { return _nodes[idx]; }
    // ------------------------------------------------------------
    inline NodeClass* node(HEClass *he) { return node(he->nodeIdx); }
    // ------------------------------------------------------------
    inline HEClass* incidentHalfEdge(const IT &idx) { return node(idx)->halfEdge; }
    // ------------------------------------------------------------
    inline void newNode(const IT &idx, NodeClass* n) { _nodes.insert(std::make_pair(idx,n)); }
    // ------------------------------------------------------------
    NodeClass* newNode(const IT &idx)
    {
      NodeClass *res = new NodeClass;
      newNode(idx,res);
      return res;
    }
    // ------------------------------------------------------------
    HEClass* newEdge(const IT &n1_idx, const IT &n2_idx)
    {
      // create halfedges
      HEClass *d12 = new HEClass(n1_idx);
      HEClass *d21 = new HEClass(n2_idx,d12);
      d12->opposite = d21;
      // insert them into array of halfedges
      _hes.push_back(d12);
      _hes.push_back(d21);
      // update orbit of n1
      NodeClass *n = node(n1_idx);
      HEClass *dn = n->halfEdge;
      n->halfEdge = d12;
      if (dn)
      {
	d12->next = dn->next;
	dn->next = d12;
      }
      // update orbit of n2
      n = node(n2_idx);
      dn = n->halfEdge;
      n->halfEdge = d21;
      if (dn)
      {
	d21->next = dn->next;
	dn->next = d21;
      }
      return d12;
    }
    // ------------------------------------------------------------
    HEClass* isNeighbor(const IT &n1_idx, const IT &n2_idx)
    {
      HEClass *he = incidentHalfEdge(n1_idx);
      if (he == NULL) return NULL;
      HEClass *end = he;
      do
      {
    	if (he->opposite->nodeIdx == n2_idx) return he;
    	he = he->next;
      } while (he != end);
      return NULL;
    }
    // ------------------------------------------------------------
    IT degree(NodeClass *n)
    {
      HEClass *he = n->halfEdge;
      if (he == NULL) return 0;
      HEClass *end = he;
      IT dg = 0;
      do
      {
	dg++;
	he = he->next;
      } while (he != end);
      return dg;
    }
    // ------------------------------------------------------------
  };

  // ------------------------------------------------------------
  template <typename IT, typename FT>
  class NodeDt
  {
  public:
    typedef FT ValueType;
    FT y;
    FT x;
    FT val;
    IT idxo;
    IT idxi;
    // ------------------------------------------------------------
    NodeDt()
      : y(-1), x(-1), val(0), idxo(-1), idxi(-1) { }
    // ------------------------------------------------------------
    NodeDt(const FT &yp, const FT &xp, const FT &value, const IT &idxout, const IT &idxin)
      : y(yp), x(xp), val(value), idxo(idxout), idxi(idxin) { }
    // ------------------------------------------------------------
    ~NodeDt() { }
    // ------------------------------------------------------------
  };
  
  template <typename IT, typename FT>
  class Node
  {
  public:
    typedef FT ValueType;
    IT dart;
    FT val;
    IT npt;
    bool border;
    IT degree;
    IT idxo;
    IT idxi;
    bool keep;
    // ------------------------------------------------------------
    Node()
      : dart(-1), val(0), npt(-1), border(false), degree(0),
	idxo(-1),idxi(-1), keep(false) { }
    // ------------------------------------------------------------
    Node(const FT &value, const IT &npoint, const bool &b)
      : dart(-1), val(value), npt(npoint), border(false), degree(0),
	idxo(-1), idxi(-1), keep(false) { }
    // ------------------------------------------------------------
    Node(const FT &value, const IT &idxout, const IT &idxin)
      : dart(-1), val(value), npt(-1), border(false), degree(0),
	idxo(idxout), idxi(idxin), keep(false) { }
    // ------------------------------------------------------------
    ~Node() { }
    // ------------------------------------------------------------
  };

  // ------------------------------------------------------------
  template <typename IT, typename FT = double>
  class Dart
  {
  public:
    IT node;
    IT alpha; // opposite dart
    IT sigma; // next dart
    FT val;
    IT idxo;
    IT idxi;
    bool keep;
    Dart() : node(-1), alpha(-1), sigma(-1), val(0), idxo(-1), idxi(-1), keep(false) { }
    Dart(const IT &nidx, const FT &value, const IT &idxout, const IT &idxin)
      : node(nidx), alpha(-1), sigma(-1), val(value),
	idxo(idxout), idxi(idxin), keep(false) { }
    ~Dart() { }
  };

  // ------------------------------------------------------------
  template <typename IT, class Node>
  class IndexedCMap
  {
  protected:
    // ------------------------------------------------------------
    std::map<IT,Node*> _nodes;
    std::vector<Dart<IT>* > _darts;
    typedef typename Node::ValueType ValueType;
  public:
    typedef typename std::map<IT,Node*>::iterator NodeIt;
    typedef typename std::vector<Dart<IT>*>::iterator DartIt;
    // ------------------------------------------------------------
    IndexedCMap()
    {
      
    }
    // ------------------------------------------------------------
    ~IndexedCMap()
    {
      for (NodeIt it = _nodes.begin(), end = _nodes.end(); it != end; ++it)
      { delete it->second; it->second = NULL; }
      for (DartIt it = _darts.begin(), end = _darts.end(); it != end; ++it)
      { delete *it; *it = NULL; }
    }
    // ------------------------------------------------------------
    Node* node()
    {
      NodeIt nit = _nodes.begin();
      if (nit == _nodes.end()) return NULL;
      return nit->second;
    }
    // ------------------------------------------------------------
    Node* node(const IT &i)
    {
      NodeIt nit = _nodes.find(i);
      if (nit == _nodes.end()) return NULL;
      return nit->second;
    }
    // ------------------------------------------------------------
    Dart<IT>* dart(const IT &i) { return _darts[i]; }
    // ------------------------------------------------------------
    IT sigmam1(const IT &didx)
    {
      if (didx == -1) return didx;
      IT didxres = didx;
      Dart<IT> *d = _darts[didx];
      while (d->sigma != didx) { didxres = d->sigma; d = _darts[didxres]; }
      return didxres;
    }
    // ------------------------------------------------------------
    const IT& addNode(const IT &nidx, Node *nd)
    {
      NodeIt nit = _nodes.find(nidx);
      if (nit == _nodes.end())
      {
	_nodes.insert(std::make_pair(nidx,nd));
	return nidx;
      }
      return nit->first;
    }
    // ------------------------------------------------------------
    const IT& addNode(const IT &nidx)
    {
      NodeIt nit = _nodes.find(nidx);
      if (nit == _nodes.end())
      {
	_nodes.insert(std::make_pair(nidx,new Node));
	return nidx;
      }
      return nit->first;
    }
    // ------------------------------------------------------------
    const IT& addNode(const IT &nidx, ValueType v,
		      const IT &idxo, const IT &idxi)
    {
      NodeIt nit = _nodes.find(nidx);
      if (nit == _nodes.end())
      {
	_nodes.insert(std::make_pair(nidx,new Node(v,idxo,idxi)));
	return nidx;
      }
      if (nit->second->val > 0)
      {
	if (v < nit->second->val)
	{
	  nit->second->val = v;
	  nit->second->idxo = idxo;
	  nit->second->idxi = idxi;
	}
      }
      else
      {
	nit->second->val = v;
	nit->second->idxo = idxo;
	nit->second->idxi = idxi;
      }
      return nit->first;
    }
    // ------------------------------------------------------------
    IT isNeighbor(const IT &nidx, const IT &idx)
    {
      NodeIt it = _nodes.find(idx);
      if (it == _nodes.end()) return -1;
      IT didx = it->second->dart;
      if (didx == -1) return -1;
      Dart<IT> *d = NULL;
      IT didxs = didx, da;
      do
      {
	d = _darts[didxs];
	da = d->alpha;
	if (_darts[da]->node == nidx) return da;
	didxs = d->sigma;
      } while (didxs != didx);
      return -1;
    }
    // ------------------------------------------------------------
    IT link(const IT &idx1, const IT &idx2)
    {
      NodeIt nit = _nodes.find(idx1);
      if (nit == _nodes.end()) return -1;
      Node *n1 = nit->second;
      nit = _nodes.find(idx2);
      if (nit == _nodes.end()) return -1;
      Node *n2 = nit->second;
      Dart<IT> *d1 = new Dart<IT>, *d2 = new Dart<IT>;
      _darts.push_back(d1); _darts.push_back(d2);
      n1->degree++; n2->degree++;
      IT id1 = _darts.size()-2, id2 = _darts.size()-1;
      d1->node = idx1; d2->node = idx2;
      d1->alpha = id2; d2->alpha = id1;
      IT dm11 = sigmam1(n1->dart), dm22 = sigmam1(n2->dart);
      if (dm11 == -1) { d1->sigma = id1; n1->dart = id1; }
      else { d1->sigma = n1->dart; _darts[dm11]->sigma = id1; }
      if (dm22 == -1) { d2->sigma = id2; n2->dart = id2; }
      else { d2->sigma = n2->dart; _darts[dm22]->sigma = id2; }
      return id1;
    }
    // ------------------------------------------------------------
    IT link(const IT &idx1, const IT &idx2, ValueType v,
	      const IT &idxo, const IT &idxi)
    {
      NodeIt nit = _nodes.find(idx1);
      if (nit == _nodes.end()) return -1;
      Node *n1 = nit->second;
      nit = _nodes.find(idx2);
      if (nit == _nodes.end()) return -1;
      Node *n2 = nit->second;
      Dart<IT> *d1 = new Dart<IT>(idx1,v,idxo,idxi), 
	*d2 = new Dart<IT>(idx2,v,idxo,idxi);
      _darts.push_back(d1); _darts.push_back(d2);
      n1->degree++; n2->degree++;
      IT id1 = _darts.size()-2, id2 = _darts.size()-1;
      d1->alpha = id2; d2->alpha = id1;
      IT dm11 = sigmam1(n1->dart), dm22 = sigmam1(n2->dart);
      if (dm11 == -1) { d1->sigma = id1; n1->dart = id1; }
      else { d1->sigma = n1->dart; _darts[dm11]->sigma = id1; }
      if (dm22 == -1) { d2->sigma = id2; n2->dart = id2; }
      else { d2->sigma = n2->dart; _darts[dm22]->sigma = id2; }
      return id1;
    }
    // ------------------------------------------------------------
    void updateDart(const IT &didx, ValueType v,
		    const IT &idxo, const IT &idxi)
    {
      Dart<IT> *d = _darts[didx];
      if (d->val > 0)
      {
	if (v < d->val)
	{
	  _darts[d->alpha]->val = d->val = v;
	  _darts[d->alpha]->idxo = d->idxo = idxo;
	  _darts[d->alpha]->idxi = d->idxi = idxi;
	}
      }
      else
      {
	_darts[d->alpha]->val = d->val = v;
	_darts[d->alpha]->idxo = d->idxo = idxo;
	_darts[d->alpha]->idxi = d->idxi = idxi;
      }
    }
    // ------------------------------------------------------------
    void degree1Nodes(std::vector<IT> &c)
    {
      for (NodeIt it = _nodes.begin(), end = _nodes.end(); it != end; ++it)
      {
	if (it->second->degree == 1) c.push_back(it->first);
	if (it->second->degree == 0) std::cerr << "zero degree node\n";
	if (it->second->degree > 2)
	  std::cerr << it->second->degree << " degree node\n";
      }
    }
    // ------------------------------------------------------------
    template <class NDt, class GridClass>
    void connectedComponents(std::vector<std::vector<NDt> > &v, const GridClass &G)
    {
      std::vector<IT> c;
      degree1Nodes(c);
      std::cerr << c.size() << " nodes of degree 1\n";
      ValueType y, x, smod, diff;
      IT d, ad, n, n1, r1, c1, rp, cp;
      
      if (c.size() > 0)
      {
	// for each connected component
	for (IT i = 0; i < c.size(); i++)
	{
	  if (_nodes[c[i]]->border) continue; // already treated
	  _nodes[c[i]]->border = true;
	  smod = 0;
	  std::vector<NDt> sv;
	  // find first edge
	  n1 = c[i];
	  d = _nodes[n1]->dart;
	  if (_nodes[n1]->keep)
	  {
	    G.point(n1,rp,cp);
	    y = rp + 0.5;
	    x = cp + 0.5;
	    NDt nk(y,x,_nodes[n1]->val,_nodes[n1]->idxo,_nodes[n1]->idxi);
	    sv.push_back(nk);
	  }
	  else
	  {
	    ad = _darts[d]->alpha;
	    n = _darts[ad]->node;
	    if (_darts[ad]->keep)
	    {
	      G.point(n1,r1,c1);
	      G.point(n,rp,cp);
	      y = (r1+rp+1)/2.0;
	      x = (c1+cp+1)/2.0;
	      NDt nk(y,x,_darts[ad]->val,_darts[ad]->idxo,_darts[ad]->idxi);
	      sv.push_back(nk);
	    }
	    else 
	      if (_nodes[n]->keep)
	      {
		G.point(n,rp,cp);
		y = rp + 0.5;
		x = cp + 0.5;
		NDt nk(y,x,_nodes[n]->val,_nodes[n]->idxo,_nodes[n]->idxi);
		sv.push_back(nk);
	      }
	      else { std::cerr << "keep error\n"; return; }
	    n1 = n;
	    d = _darts[ad]->sigma;
	  }
	  // continue
	  do
	  {
	    ad = _darts[d]->alpha;
	    n = _darts[ad]->node;
	    _nodes[n]->border = true;
	    if (_darts[ad]->keep)
	    {
	      G.point(n1,r1,c1);
	      G.point(n,rp,cp);
	      y = (r1+rp+1)/2.0;
	      x = (c1+cp+1)/2.0;
	      diff = _darts[ad]->val - sv[sv.size()-1].val;
	      if (smod == 0) smod = diff;
	      if (smod > 0 && diff < 0)
	      {
	      	std::cerr << "max " << sv[sv.size()-1].y+1 << " " << sv[sv.size()-1].x+1 << " " <<sv[sv.size()-1].val<< "\n";
	      	smod = diff;
	      }
	      if (smod < 0 && diff > 0)
	      {
	      	std::cerr << "min " << sv[sv.size()-1].y+1 << " " << sv[sv.size()-1].x+1 << " " <<sv[sv.size()-1].val<< "\n";
	      	smod = diff;
	      }
	      NDt nk(y,x,_darts[ad]->val,_darts[ad]->idxo,_darts[ad]->idxi);
	      sv.push_back(nk);
	    }
	    else if (_nodes[n]->keep)
	      {
		G.point(n,rp,cp);
		y = rp + 0.5;
		x = cp + 0.5;
		diff = _nodes[n]->val - sv[sv.size()-1].val;
		if (smod == 0) smod = diff;
		if (smod > 0 && diff < 0)
		{
		  std::cerr << "max " << sv[sv.size()-1].y+1 << " " << sv[sv.size()-1].x+1 << " " <<sv[sv.size()-1].val<< "\n";
		  smod = diff;
		}
		if (smod < 0 && diff > 0)
		{
		  std::cerr << "min " << sv[sv.size()-1].y+1 << " " << sv[sv.size()-1].x+1 << " " <<sv[sv.size()-1].val<< "\n";
		  smod = diff;
		}
		NDt nk(y,x,_nodes[n]->val,_nodes[n]->idxo,_nodes[n]->idxi);
		sv.push_back(nk);
	      }
	    d = _darts[ad]->sigma;
	    n1 = n;
	  } while (d != ad);
	  v.push_back(sv);
	}
      }
    }
    // ------------------------------------------------------------
    
  };
}

#endif
