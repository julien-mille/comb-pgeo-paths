// =========================================================================
/* Copyright 2011 Sebastien Bougleux

   This file is part of Geodesics4Meshes.

   Geodesics4Meshes is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   Geodesics4Meshes is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU General Public License,
   and a copy of the GNU Lesser General Public License, along with
   GeoComp. If not, see <http://www.gnu.org/licenses/>.
*/
// =========================================================================
/**
 * @file minheap.hh
 * @authors Sebastien Bougleux, Julien Mille
 * @date May 2015
 */
// =========================================================================
#ifndef __MIN_HEAP_HH__
#define __MIN_HEAP_HH__
#include <iostream>
#include <stdlib.h>
#include "utility.hh"

namespace GeoComp
{
  // ================================================================
  /**
   * @class MinHeap
   * @brief Min heap with updates (using a binary tree and arrays of fixed size).
   * @param Container Type of data container (double* by default).
   * @param IT Signed integer type for array indices (long by default).
   * @authors Sebastien Bougleux
   * @date December 2011
   */
  template<class Container = double, class IT = long>
  class MinHeap
  {
  protected:
    /**
     * @brief Fixed (maximum) number of elements in the heap.
     */
    IT _size;
    /**
     * @brief Number of elements in the heap.
     */
    IT _cur_size;
    /**
     * @brief Tree structure (an array of size _size) such that _index[i] provides the index of the element of _data[i] in the sorted list _btree.
     */
    IT *_index;
    /**
     * @brief List of the element indices in _data (part of the binary tree structure) sorted by their values in _data.
     */
    IT *_btree;
    /**
     * @brief The data container, external to the class. It must support operator [] to acces the data, and operators < and > to compare any data pair.
     */
    Container *_data;

    //================================================================
    /**
     * @brief Initialize the heap.
     */
    inline void initialize()
    {
      for (IT i = 0; i < _size; i++)
        _index[i] = -1;
    }
    //================================================================
    /**
     * @brief Initialize the heap.
     */
    inline void reset()
    {
      IT j;
      for (j=0; j<_cur_size; j++)
        _index[_btree[j]] = -1;
      _cur_size = 0;
    }
    //================================================================
  public:
    //================================================================
    /**
     * @brief Constructor.
     * @param[in] s Size of arrays.
     * @param[in] u Data container.
     */
    MinHeap(IT s, Container *c)
      : _size(s), _cur_size(0), _index(new IT[s]), _btree(new IT[s]), _data(c)
    { initialize(); }
    //================================================================
    /**
     * @brief Destructor.
     */
    ~MinHeap() { delete[] _index; delete[] _btree; }
    //================================================================
    /**
     *  @brief Test if the heap is empty.
     *  @return true if it is empty, or false else.
     */
    inline bool empty() { return (_cur_size == 0); }
    //================================================================
    /**
     * @brief Clear the min heap.
     */
    inline void clear()
    {
      if (!empty()) reset();
    }
    //================================================================
    /**
     * @brief Set the data array.
     * @param[in] d Data array.
     */
    inline void setData(Container *d) { _data = d; }
    //================================================================
    /**
     * @brief Push a new element in the heap.
     * @param[in] idxDataInArray Index of the element in the data array.
     */
    IT push(IT idxInDataArray)
    {
      if (_cur_size==0)
      {
        _index[idxInDataArray] = 0;
        _btree[0] = idxInDataArray;
      }
      else {
        IT idx = _cur_size, parent_idx=(idx-1)/2;
        IT idx_temp;
        _index[idxInDataArray] = idx;
        _btree[idx] = idxInDataArray;

        while (idx>0 && _data[_btree[parent_idx]]>_data[_btree[idx]])
        {
            _index[_btree[parent_idx]] = idx;
            _index[_btree[idx]] = parent_idx;

            idx_temp = _btree[parent_idx];
            _btree[parent_idx] = _btree[idx];
            _btree[idx] = idx_temp;
            idx = parent_idx;
            parent_idx = (parent_idx-1)/2;
        }
      }
      _cur_size++;

      return _cur_size;
    }
    //================================================================
    /**
     * @brief First element of the heap (the minimum element).
     */
    IT top() {
        assert(_cur_size>0);
        return *_btree;
    }
    //================================================================
    /**
     * @brief Pop the first element of the heap (the minimum element).
     */
    IT pop()
    {
      if (_cur_size<=0)
      {
        std::cout<<"pop() : heap is empty\n";
        exit(-1);
      }

      IT first = *_btree;
      IT idx, leftchild_idx, rightchild_idx, idx_temp;
        bool permutation = true;

       _index[_btree[0]] = -1;
       _index[_btree[_cur_size-1]] = 0;

       _btree[0] = _btree[_cur_size-1];
       _cur_size--;

        idx = 0;
        leftchild_idx = 1;
        rightchild_idx = 2;

        while (leftchild_idx<_cur_size && permutation==true)
        {
            // std::cout<<_btree[idx]<<"\n";

            if (!(_btree[idx]>=0 && _btree[idx]<_size))
                std::cout<<"idx wrong\n";
            if (!(_btree[leftchild_idx]>=0 && _btree[leftchild_idx]<_size))
                std::cout<<"leftchild_idx wrong\n";
            if (!(_btree[rightchild_idx]>=0 && _btree[rightchild_idx]<_size))
                std::cout<<"rightchild_idx wrong\n";

            if (rightchild_idx<_cur_size && _data[_btree[rightchild_idx]]<_data[_btree[idx]] &&
                _data[_btree[rightchild_idx]]<_data[_btree[leftchild_idx]])
            {
                _index[_btree[rightchild_idx]] = idx;
                _index[_btree[idx]] = rightchild_idx;

                idx_temp = _btree[rightchild_idx];
                _btree[rightchild_idx] = _btree[idx];
                _btree[idx] = idx_temp;

                idx = rightchild_idx;
            }
            else if (_data[_btree[leftchild_idx]]<_data[_btree[idx]])
            {
                _index[_btree[leftchild_idx]] = idx;
                _index[_btree[idx]] = leftchild_idx;

                idx_temp = _btree[leftchild_idx];
                _btree[leftchild_idx] = _btree[idx];
                _btree[idx] = idx_temp;

                idx = leftchild_idx;
            }
            else
                permutation = false;

            leftchild_idx = 2*idx+1;
            rightchild_idx = 2*idx+2;
        }
      return first;
    }
    //================================================================
    /**
     * @brief Update the position, in the heap, of a given element.
     * @param[in] idxInDataArray Index of the element in the data array.
     */
    void update(IT idxInDataArray)
    {
      IT idx = _index[idxInDataArray], parent_idx=(idx-1)/2;
      IT idx_temp;
        while (idx>0 && _data[_btree[parent_idx]]>_data[_btree[idx]])
        {
            _index[_btree[parent_idx]] = idx;
            _index[_btree[idx]] = parent_idx;

            idx_temp = _btree[parent_idx];
            _btree[parent_idx] = _btree[idx];
            _btree[idx] = idx_temp;

            idx = parent_idx;
            parent_idx = (idx-1)/2;
        }
    }
    //================================================================
    /**
     * @brief Update the position, in the heap, of a given element.
     * @param[in] idxInDataArray Index of the element in the data array.
     */
    void printData()
    {
      IT i;
      std::cout << "index" << std::endl;
      for (i = 0; i < _size; i++) std::cout << _index[i] << " ";
      std::cout << std::endl << "btree" << std::endl;
      for (i = 0; i < _size; i++) std::cout << _btree[i] << " ";
      std::cout << std::endl << std::endl;
    }
  };
}

#endif
