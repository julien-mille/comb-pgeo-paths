/*
Copyright 2010 Julien Mille

This file is part of CombinationOfPiecewiseGeodesicPaths.

CombinationOfPiecewiseGeodesicPaths is free software: you can redistribute
it and/or modify it under the terms of the GNU Lesser General Public License
as published by the Free Software Foundation, either version 3 of the License,
or (at your option) any later version.

CombinationOfPiecewiseGeodesicPaths is distributed in the hope that it will
be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser
General Public License for more details.

You should have received a copy of the GNU General Public License,
and a copy of the GNU Lesser General Public License, along with
CombinationOfPiecewiseGeodesicPaths. If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef _ARRAYND_H_
#define _ARRAYND_H_

#include <limits>
#include <iostream>
#include "couple.h"
#include "triplet.h"

// #define ARRAYND_RUNTIME_CHECK

using namespace std;

template <class T> class CArray1D;
template <class T> class CArray2D;
template <class T> class CArray3D;

template <class T> class CArray1D
{
  // Static members
#ifdef ARRAYND_RUNTIME_CHECK
  protected:
	// Element returned in case of invalid access with operator [] or function GetElementInterpolate()
	static T elementError;
#endif

  // Members
  protected:
	int iSize;
	T *pElements;

  // Member functions
  public:
	inline CArray1D();
	inline CArray1D(int);
	inline CArray1D(const CArray1D<T> &);

	inline ~CArray1D();

	inline int GetSize() const {return iSize;}
	inline const T *GetBuffer() const {return pElements;}
	inline T *GetBuffer() {return pElements;}

	inline bool Init(int);
	inline bool Init(int, const T *);
	template <class S> bool InitCast(const CArray1D<S> &);
	inline void Empty();
	inline void Fill(const T &);

	inline CArray1D &operator =(const CArray1D<T> &);

	// Access single element (read/write)
	inline T &operator [](int index) const
	{
		#ifdef ARRAYND_RUNTIME_CHECK
		if (index<0 || index>=iSize)
		{
			cerr<<"ERROR in CArray1D<T>::operator [](int): accessing element at index "<<index<<" out of range [0.."<<iSize-1<<"]"<<endl;
			return elementError;
		}
		#endif
		return pElements[index];
	}

	// Linear interpolation between elements (read only)
	inline T GetElementInterpolate(float fIndex) const
	{
		float fDiff;
		int iIndex;
		T elementInterpolate;

		#ifdef ARRAYND_RUNTIME_CHECK
		if (fIndex<0.0f || fIndex>=(float)(iSize-1))
		{
			cerr<<"ERROR in CArray1D<T>::GetElementInterpolate(float): accessing element at index "<<fIndex<<" out of range [0.."<<iSize-1<<"["<<endl;
			return elementError;
		}
		#endif

		iIndex = (int)fIndex;
		fDiff = fIndex-floor(fIndex);
		elementInterpolate = (T)((1.0f-fDiff)*pElements[iIndex] + fDiff*pElements[iIndex+1]);

		return elementInterpolate;
	}

	CArray1D<T> &operator +=(const CArray1D<T> &);
	CArray1D<T> &operator -=(const CArray1D<T> &);
	template <class S> CArray1D<T> &operator *=(const S &);
	template <class S> CArray1D<T> &operator /=(const S &);

	bool operator ==(const CArray1D<T> &) const;
	CArray1D<T> operator +(const CArray1D<T> &) const;
	CArray1D<T> operator -() const;
	CArray1D<T> operator -(const CArray1D<T> &) const;
	template <class S> CArray1D<T> operator *(const S &) const;

	template <class U, class S> friend CArray1D<U> operator *(const S &, const CArray1D<U> &);
	template <class S> CArray1D<T> operator /(const S &) const;

	T Sum() const;
};

#ifdef ARRAYND_RUNTIME_CHECK
// Element returned in case of invalid access with operator [] or function GetElementInterpolate()
template <class T> T CArray1D<T>::elementError;
#endif

template <class T> CArray1D<T>::CArray1D()
{
	iSize = 0;
	pElements = NULL;
}

template <class T> CArray1D<T>::CArray1D(int size)
{
	iSize = size;
	pElements = new T[iSize];
}

template <class T> CArray1D<T>::CArray1D(const CArray1D<T> &array)
{
	T *pElementsTemp;
	const T *pElementsTemp2;
	int i;

	iSize = array.iSize;
	pElements = new T[iSize];

	pElementsTemp = pElements;
	pElementsTemp2 = array.pElements;
	for (i=0;i<iSize;i++)
		*pElementsTemp++ = *pElementsTemp2++;
}

template <class T> CArray1D<T>::~CArray1D()
{
	Empty();
}

template <class T> bool CArray1D<T>::Init(int size)
{
	if (pElements!=NULL)
		delete[] pElements;

	iSize = size;

	pElements = new T[iSize];
	if (pElements!=NULL)
		return true;
	else
		return false;
}

template <class T> bool CArray1D<T>::Init(int size, const T *ptr)
{
	if (pElements!=NULL)
		delete[] pElements;

	iSize = size;

	pElements = new T[iSize];
	if (pElements!=NULL)
	{
		T *pElementsTemp = pElements;
		const T *ptrTemp = ptr;
		int i;

		for (i=0;i<iSize;i++)
			*pElementsTemp++ = *ptrTemp++;
		return true;
	}
	else
		return false;
}

template <class T> template <class S> bool CArray1D<T>::InitCast(const CArray1D<S> &array)
{
	if (pElements!=NULL)
		delete[] pElements;

	iSize = array.GetSize();

	pElements = new T[iSize];
	if (pElements!=NULL)
	{
		T *pElementsThis = pElements;
		const S *pElementsParam = array.GetBuffer();
		int i;

		for (i=0;i<iSize;i++)
		{
			*pElementsThis = (T)(*pElementsParam);
			pElementsThis++;
			pElementsParam++;
		}
		return true;
	}
	else
		return false;
}

template <class T> void CArray1D<T>::Empty()
{
	if (pElements!=NULL)
		delete[] pElements;

	iSize = 0;
	pElements = NULL;
}

template <class T> void CArray1D<T>::Fill(const T &t)
{
	int i;
	T *pElementsTemp = pElements;

	for (i=0;i<iSize;i++)
	{
		*pElementsTemp = t;
		pElementsTemp++;
	}
}

template <class T> CArray1D<T> &CArray1D<T>::operator =(const CArray1D<T> &array)
{
	int i;

	// Reallocate if necessary
	if (iSize!=array.iSize)
	{
		if (pElements!=NULL)
		{
			delete[] pElements;
			pElements = NULL;
		}
		iSize = array.iSize;
		pElements = new T[iSize];
	}

	// Copy elements
	if (iSize!=0)
	{
		T *pElementsTemp;
		const T *pElementsTemp2;
		pElementsTemp = pElements;
		pElementsTemp2 = array.pElements;

		for (i=0;i<iSize;i++)
			*pElementsTemp++ = *pElementsTemp2++;
	}
	return *this;
}

template <class T> CArray1D<T> &CArray1D<T>::operator +=(const CArray1D<T> &array)
{
	if (iSize==array.iSize)
	{
		T *pElementsTemp;
		const T *pElementsTemp2;
		int i;

		pElementsTemp = pElements;
		pElementsTemp2 = array.pElements;

		for (i=0;i<iSize;i++)
			*pElementsTemp++ += *pElementsTemp2++;
		return *this;
	}
	else throw 0;
}

template <class T> CArray1D<T> &CArray1D<T>::operator -=(const CArray1D<T> &array)
{
	if (iSize==array.iSize)
	{
		T *pElementsTemp;
		const T *pElementsTemp2;
		int i;

		pElementsTemp = pElements;
		pElementsTemp2 = array.pElements;

		for (i=0;i<iSize;i++)
			*pElementsTemp++ -= *pElementsTemp2++;
		return *this;
	}
	else throw 0;
}

template <class T> template <class S> CArray1D<T> &CArray1D<T>::operator *=(const S &t)
{
	T *pElementsTemp = pElements;
	int i;

	for (i=0;i<iSize;i++)
		*pElementsTemp++ *= t;
	return *this;
}

template <class T> template <class S> CArray1D<T> &CArray1D<T>::operator /=(const S &t)
{
	T *pElementsTemp = pElements;
	int i;

	for (i=0;i<iSize;i++)
		*pElementsTemp++ /= t;
	return *this;
}

template <class T> bool CArray1D<T>::operator ==(const CArray1D<T> &array) const
{
	if (iSize!=array.iSize)
		return false;

	const T *pElementsTemp, *pElementsTemp2;
	int i;

	pElementsTemp = this->pElements;
	pElementsTemp2 = array.pElements;
	for (i=0;i<this->iSize && *pElementsTemp==*pElementsTemp2; i++)
	{
		pElementsTemp++;
		pElementsTemp2++;
	}
	if (i==iSize) return true;
	else return false;
}

template <class T> CArray1D<T> CArray1D<T>::operator +(const CArray1D<T> &array) const
{
	if (iSize==array.iSize)
	{
		CArray1D<T> arrayResult(iSize);
		const T *pElementsTemp, *pElementsTemp2;
		T *pElementsTempResult;
		int i;

		pElementsTemp = pElements;
		pElementsTemp2 = array.pElements;
		pElementsTempResult = arrayResult.pElements;

		for (i=0;i<iSize;i++)
		{
			*pElementsTempResult = *pElementsTemp + *pElementsTemp2;
			pElementsTemp++;
			pElementsTemp2++;
			pElementsTempResult++;
		}
		return arrayResult;
	}
	else throw 0;
}

template <class T> CArray1D<T> CArray1D<T>::operator -() const
{
	CArray1D<T> arrayResult(iSize);
	const T *pElementsTemp;
	T *pElementsTempResult;
	int i;

	pElementsTemp = pElements;
	pElementsTempResult = arrayResult.pElements;

	for (i=0;i<iSize;i++)
	{
		*pElementsTempResult = -(*pElementsTemp);
		pElementsTemp++;
		pElementsTempResult++;
	}
	return arrayResult;
}

template <class T> CArray1D<T> CArray1D<T>::operator -(const CArray1D<T> &array) const
{
	if (iSize==array.iSize)
	{
		CArray1D<T> arrayResult(iSize);
		const T *pElementsTemp, *pElementsTemp2;
		T *pElementsTempResult;
		int i;

		pElementsTemp = pElements;
		pElementsTemp2 = array.pElements;
		pElementsTempResult = arrayResult.pElements;

		for (i=0;i<iSize;i++)
		{
			*pElementsTempResult = (*pElementsTemp)-(*pElementsTemp2);
			pElementsTemp++;
			pElementsTemp2++;
			pElementsTempResult++;
		}
		return arrayResult;
	}
	else throw 0;
}

template <class T> template <class S> CArray1D<T> CArray1D<T>::operator *(const S &t) const
{
	CArray1D<T> arrayResult(iSize);
	const T *pElementsTemp;
	T *pElementsTempResult;
	int i;

	pElementsTemp = pElements;
	pElementsTempResult = arrayResult.pElements;
	for (i=0;i<iSize;i++)
	{
		*pElementsTempResult = (*pElementsTemp)*t;
		pElementsTemp++;
		pElementsTempResult++;
	}
	return arrayResult;
}

template <class T, class S> CArray1D<T> operator *(const S &t, const CArray1D<T> &v)
{
	CArray1D<T> arrayResult(v.iSize);
	const T *pElementsTemp;
	T *pElementsTempResult;
	int i;

	pElementsTemp = v.pElements;
	pElementsTempResult = arrayResult.pElements;
	for (i=0;i<v.iSize;i++)
	{
		*pElementsTempResult = (*pElementsTemp)*t;
		pElementsTemp++;
		pElementsTempResult++;
	}
	return arrayResult;
}

template <class T> template <class S> CArray1D<T> CArray1D<T>::operator /(const S &t) const
{
	CArray1D<T> arrayResult(iSize);
	const T *pElementsTemp;
	T *pElementsTempResult;
	int i;

	pElementsTemp = pElements;
	pElementsTempResult = arrayResult.pElements;
	for (i=0;i<iSize;i++)
	{
		*pElementsTempResult = (*pElementsTemp)/t;
		pElementsTemp++;
		pElementsTempResult++;
	}
	return arrayResult;
}

template <class T> T CArray1D<T>::Sum() const
{
    const T *pElementsTemp;
	T sum = (T)0;
	int i;

	pElementsTemp = pElements;
	for (i=0;i<iSize;i++)
	{
	    sum += *pElementsTemp;
		pElementsTemp++;
	}
	return sum;
}

template <class T> class CArray2DIterator
{
  friend class CArray2D<T>;

  // Attributes
  protected:
	T *pCurrentElementY, *pCurrentElement;
	CCouple<int> ptCurrent, ptMin, ptMax;
	CArray2D<T> *pArrayParent;

  public:
	inline CArray2DIterator()
	{
		pArrayParent = NULL;
		pCurrentElementY = pCurrentElement = NULL;
	}
	T *operator ++() // Prefix operator ++
	{
		T *pRet;

		ptCurrent.x++;
		pCurrentElement++;
		if (ptCurrent.x>ptMax.x)
		{
			ptCurrent.x = ptMin.x;
			ptCurrent.y++;

			pCurrentElementY += pArrayParent;
			pCurrentElement = pCurrentElementY + ptMin.x;
		}
		pRet = pCurrentElement;
		return pRet;
	}
	T *operator ++(int) // Postfix operator ++
	{
		T *pRet;

		pRet = pCurrentElement;

		ptCurrent.x++;
		pCurrentElement++;
		if (ptCurrent.x>ptMax.x)
		{
			ptCurrent.x = ptMin.x;
			ptCurrent.y++;

			pCurrentElementY += pArrayParent->iWidth;
			pCurrentElement = pCurrentElementY + ptMin.x;
		}
		return pRet;
	}

	inline CCouple<int> GetPosition() {return ptCurrent;}
	inline T &Element() {return *pCurrentElement;}
	inline T *ElementPtr() {return pCurrentElement;}
	inline operator T *() {return pCurrentElement;}
	inline int End() {return (ptCurrent.y>ptMax.y);}
};

template <class T> class CArray2D : protected CArray1D<T>
{
	friend class CArray2DIterator<T>;

  protected:
	int iWidth, iHeight;

  // Member functions
  public:
	inline CArray2D();
	inline CArray2D(int, int);
	inline CArray2D(const CCouple<int> &);
	inline CArray2D(const CArray2D<T> &);

	inline ~CArray2D();

	inline int GetWidth() const {return iWidth;}
	inline int GetHeight() const {return iHeight;}
	inline CCouple<int> GetSize() const {return CCouple<int>(iWidth, iHeight);}
	inline int GetOffset(int x, int y) const {return y*iWidth+x;}
	inline int GetOffset(const CCouple<int> &coord) const {return coord.y*iWidth + coord.x;}

	inline const T *GetBuffer() const {return this->pElements;}
    inline T *GetBuffer() {return this->pElements;}

	bool Init(int, int);
	bool Init(const CCouple<int> &);
	template <class S> bool InitCast(const CArray2D<S> &);
	void Empty();
	void Fill(const T &);
	void FillBoundary(const T &);

	CArray2D<T> &operator =(const CArray2D<T> &);

	// Access single element (read/write)
	inline T &Element(int x, int y) const
	{
		#ifdef ARRAYND_RUNTIME_CHECK
		if (x<0 || x>=iWidth || y<0 || y>=iHeight)
		{
			cerr<<"ERROR in CArray2D<T>::Element(int, int): accessing element ("<<x<<","<<y<<") out of range [0.."<<iWidth-1<<"]x[0.."<<iHeight-1<<"]"<<endl;
			return CArray1D<T>::elementError;
		}
		#endif
		return this->pElements[y*iWidth+x];
	}
	inline T &Element(const CCouple<int> &coord) const
	{
		#ifdef ARRAYND_RUNTIME_CHECK
		if (coord.x<0 || coord.x>=iWidth || coord.y<0 || coord.y>=iHeight)
		{
			cerr<<"ERROR in CArray2D<T>::Element(const CCouple<int> &): accessing element ("<<coord.x<<","<<coord.y<<") out of range [0.."<<iWidth-1<<"]x[0.."<<iHeight-1<<"]"<<endl;
			return CArray1D<T>::elementError;
		}
		#endif
		return this->pElements[coord.y*iWidth+coord.x];
	}

	// Linear interpolation between elements (read only)
	inline T GetElementInterpolate(float x, float y) const
	{
		float dx, dy;
		int xi, yi;
		const T *pElementTemp;
		T elementInterpolate;

		#ifdef ARRAYND_RUNTIME_CHECK
		if (x<0.0f || x>=(float)(iWidth-1) || y<0.0f || y>=(float)(iHeight-1))
		{
			cerr<<"ERROR in CArray2D<T>::GetElementInterpolate(float,float): accessing pixel ("<<x<<","<<y<<") out of range [0.."<<iWidth-1<<"[x[0.."<<iHeight-1<<"["<<endl;
			return elementInterpolate;
		}
		#endif

		xi = (int)x;
		yi = (int)y;
		dx = x-floor(x);
		dy = y-floor(y);

		// Get address of nearest element with lower integer coordinates
		pElementTemp = this->pElements + yi*this->iWidth + xi;

		elementInterpolate = (T)(
			(1.0f-dx)*(1.0f-dy) * pElementTemp[0] +
			dx*(1.0f-dy)        * pElementTemp[1] +
			(1.0f-dx)*dy        * pElementTemp[iWidth] +
			dx*dy               * pElementTemp[iWidth+1]
		);
		return elementInterpolate;
	}
	inline T GetElementInterpolate(const CCouple<float> &coord) const {return GetElementInterpolate(coord.x,coord.y);}

	// Iterator
	inline CArray2DIterator<T> GetIterator()
	{
		CArray2DIterator<T> iterator;

		iterator.pArrayParent = this;
		iterator.pCurrentElementY = iterator.pCurrentElement = this->pElements;
		iterator.ptCurrent.Set(0,0);
		iterator.ptMin.Set(0,0);
		iterator.ptMax.Set(iWidth-1,iHeight-1);

		return iterator;
	}
	inline CArray2DIterator<T> GetIterator(const CCouple<int> &ptMin, const CCouple<int> &ptMax)
	{
		CArray2DIterator<T> iterator;
		CCouple<int> ptLimitLow(0), ptLimitHigh(GetSize()-CCouple<int>(1));

		iterator.pArrayParent = this;
		iterator.pCurrentElementY = this->pElements + ptMin.y*iWidth;
		iterator.pCurrentElement = iterator.pCurrentElementY + ptMin.x;
		iterator.ptCurrent = iterator.ptMin = ptMin;
		iterator.ptMax = ptMax;

		if (ptMax.IsInRange(ptLimitLow, ptLimitHigh)
			&& ptMin.IsInRange(ptLimitLow, ptMax))
		{
			iterator.ptCurrent = iterator.ptMin = ptMin;
			iterator.ptMax = ptMax;
		}
		else {
			iterator.ptMax.y = 0;
			iterator.ptCurrent.y = 1;
			cerr<<"ERROR: CArray2D<T>::GetIterator()"<<endl;
		}

		return iterator;
	}

	CArray2D<T> &operator +=(const CArray2D<T> &);
	CArray2D<T> &operator -=(const CArray2D<T> &);
	template <class S> CArray2D<T> &operator *=(const S &);
	template <class S> CArray2D<T> &operator /=(const S &);

	bool operator ==(const CArray2D<T> &) const;
	CArray2D<T> operator +(const CArray2D<T> &) const;
	CArray2D<T> operator -() const;
	CArray2D<T> operator -(const CArray2D<T> &) const;
	template <class S> CArray2D<T> operator *(const S &) const;
	template <class U, class S> friend CArray2D<U> operator *(const S &, const CArray2D<U> &);
	template <class S> CArray2D<T> operator /(const S &) const;

    T Sum() const;

	CArray1D<int> GetNeighborOffsetsConnex4() const;
	CArray1D<int> GetNeighborOffsetsConnex8() const;

	// Convolution
	template <class S> CArray2D<T> Convolve(const CArray2D<S> &) const;
	template <class S> CArray2D<T> Convolve(const CArray1D<CCouple<int> > &, const CArray1D<S> &) const;

    template <class S> CArray2D<T> ConvolveMarked(const CArray2D<S> &, const CArray2D<bool> &) const;
	template <class S> CArray2D<T> ConvolveMarked(const CArray1D<CCouple<int> > &, const CArray1D<S> &, const CArray2D<bool> &) const;

    // Sub array
    CArray2D<T> SubArray(const CCouple<int> &, const CCouple<int> &) const;
};

template <class T> CArray2D<T>::CArray2D():CArray1D<T>()
{
	iWidth = iHeight = 0;
}

template <class T> CArray2D<T>::CArray2D(int width, int height):CArray1D<T>(width*height)
{
	iWidth = width;
	iHeight = height;
}

template <class T> CArray2D<T>::CArray2D(const CCouple<int> &iDimensions):CArray1D<T>(iDimensions.x*iDimensions.y)
{
	iWidth = iDimensions.x;
	iHeight = iDimensions.y;
}

template <class T> CArray2D<T>::CArray2D(const CArray2D<T> &array):CArray1D<T>(array)
{
	iWidth = array.iWidth;
	iHeight = array.iHeight;
}

template <class T> CArray2D<T>::~CArray2D()
{
	Empty();
}

template <class T> bool CArray2D<T>::Init(int width, int height)
{
	iWidth = width;
	iHeight = height;
	return CArray1D<T>::Init(iWidth*iHeight);
}

template <class T> bool CArray2D<T>::Init(const CCouple<int> &iDimensions)
{
	iWidth = iDimensions.x;
	iHeight = iDimensions.y;
	return CArray1D<T>::Init(iWidth*iHeight);
}

template <class T> template <class S> bool CArray2D<T>::InitCast(const CArray2D<S> &array)
{
	iWidth = array.GetWidth();
	iHeight = array.GetHeight();

	if (this->pElements!=NULL)
		delete[] this->pElements;

	this->iSize = iWidth*iHeight;

	this->pElements = new T[this->iSize];
	if (this->pElements!=NULL)
	{
		T *pElementsThis = this->pElements;
		const S *pElementsParam = array.GetBuffer();
		int i;

		for (i=0;i<this->iSize;i++)
		{
			*pElementsThis = (T)(*pElementsParam);
			pElementsThis++;
			pElementsParam++;
		}
		return true;
	}
	else
		return false;
}

template <class T> void CArray2D<T>::Empty()
{
	iWidth = iHeight = 0;
	CArray1D<T>::Empty();
}

template <class T> void CArray2D<T>::Fill(const T &t)
{
	CArray1D<T>::Fill(t);
}

template <class T> void CArray2D<T>::FillBoundary(const T &t)
{
	if (this->pElements==NULL)
        return;

    T *pElementsTemp = this->pElements;
    int x, y;

    // Top row
    pElementsTemp = this->pElements;
    for (x=0;x<iWidth;x++)
        *pElementsTemp++ = t;

    // Bottom row
    pElementsTemp = this->pElements + (iHeight-1)*iWidth;
    for (x=0;x<iWidth;x++)
        *pElementsTemp++ = t;

    // Left column
    pElementsTemp = this->pElements;
    for (y=0;y<iHeight;y++)
    {
        *pElementsTemp = t;
        pElementsTemp+=iWidth;
    }

    // Right column
    pElementsTemp = this->pElements + iWidth-1;
    for (y=0;y<iHeight;y++)
    {
        *pElementsTemp = t;
        pElementsTemp+=iWidth;
    }
}

template <class T> CArray2D<T> &CArray2D<T>::operator =(const CArray2D<T> &array)
{
	iWidth = array.iWidth;
	iHeight = array.iHeight;

	CArray1D<T>::operator =(array);

	return *this;
}

template <class T> CArray2D<T> &CArray2D<T>::operator +=(const CArray2D<T> &array)
{
	if (iWidth==array.iWidth && iHeight==array.iHeight)
	{
		CArray1D<T>::operator +=(array);
		return *this;
	}
	else
		throw 0;
}

template <class T> CArray2D<T> &CArray2D<T>::operator -=(const CArray2D<T> &array)
{
	if (iWidth==array.iWidth && iHeight==array.iHeight)
	{
		CArray1D<T>::operator -=(array);
		return *this;
	}
	else
		throw 0;
}

template <class T> template <class S> CArray2D<T> &CArray2D<T>::operator *=(const S &t)
{
	CArray1D<T>::operator *=(t);
	return *this;
}

template <class T> template <class S> CArray2D<T> &CArray2D<T>::operator /=(const S &t)
{
	CArray1D<T>::operator /=(t);
	return *this;
}

template <class T> bool CArray2D<T>::operator ==(const CArray2D<T> &array) const
{
	if (iWidth!=array.iWidth || iHeight!=array.iHeight)
		return false;

	const T *pElementsTemp, *pElementsTemp2;
	int i;

	pElementsTemp = this->pElements;
	pElementsTemp2 = array.pElements;
	for (i=0;i<this->iSize && *pElementsTemp==*pElementsTemp2; i++)
	{
		pElementsTemp++;
		pElementsTemp2++;
	}
	if (i==this->iSize) return true;
	else return false;
}

template <class T> CArray2D<T> CArray2D<T>::operator +(const CArray2D<T> &array) const
{
	if (iWidth==array.iWidth && iHeight==array.iHeight)
	{
		CArray2D<T> arrayRes(*this);
		arrayRes.CArray1D<T>::operator +=(array);
		return arrayRes;
	}
	else
		throw 0;
}

template <class T> CArray2D<T> CArray2D<T>::operator -() const
{
	CArray2D<T> arrayResult(iWidth, iHeight);
	const T *pElementsTemp;
	T *pElementsTempResult;
	int i;

	pElementsTemp = this->pElements;
	pElementsTempResult = arrayResult.pElements;

	for (i=0;i<this->iSize;i++)
	{
		*pElementsTempResult = -(*pElementsTemp);
		pElementsTemp++;
		pElementsTempResult++;
	}
	return arrayResult;
}

template <class T> CArray2D<T> CArray2D<T>::operator -(const CArray2D<T> &array) const
{
	if (iWidth==array.iWidth && iHeight==array.iHeight)
	{
		CArray2D<T> arrayRes(*this);
		arrayRes.CArray1D<T>::operator -=(array);
		return arrayRes;
	}
	else
		throw 0;
}

template <class T> template <class S> CArray2D<T> CArray2D<T>::operator *(const S &t) const
{
	CArray2D<T> arrayRes(*this);
	arrayRes.CArray1D<T>::operator *=(t);
	return arrayRes;
}

template <class T, class S> CArray2D<T> operator *(const S &t, const CArray2D<T> &array)
{
	CArray2D<T> arrayRes(array);
	arrayRes.CArray1D<T>::operator *=(t);
	return arrayRes;
}

template <class T> template <class S> CArray2D<T> CArray2D<T>::operator /(const S &t) const
{
	CArray2D<T> arrayRes(*this);
	arrayRes.CArray1D<T>::operator /=(t);
	return arrayRes;
}

template <class T> T CArray2D<T>::Sum() const
{
	return CArray1D<T>::Sum();
}

template <class T>  CArray1D<int> CArray2D<T>::GetNeighborOffsetsConnex4() const
{
	CArray1D<int> arrayOffsets(4);

	arrayOffsets[0] = -1;
	arrayOffsets[1] = 1;
	arrayOffsets[2] = -iWidth;
	arrayOffsets[3] = iWidth;

	return arrayOffsets;
}

template <class T>  CArray1D<int> CArray2D<T>::GetNeighborOffsetsConnex8() const
{
	CArray1D<int> arrayOffsets(8);
	CCouple<int> pi;
	int i;

	i = 0;
	for (pi.y=-1;pi.y<=1;pi.y++)
	{
		for (pi.x=-1;pi.x<=1;pi.x++)
		{
			if (pi.x!=0 || pi.y!=0)
			{
				arrayOffsets[i++] = pi.y*iWidth + pi.x;
			}
		}
	}
	return arrayOffsets;
}

// Convolution with a centered mask
template <class T> template <class S> CArray2D<T> CArray2D<T>::Convolve(const CArray2D<S> &arrayCenteredMask) const
{
	CArray1D<CCouple<int> > arrayNeighbors;
	CArray1D<S> arrayMask;
	CCouple<int> neighbor, *pNeighbor;
	S *pMask;
	const S *pCenteredMask;
	int iSizeMask;

	iSizeMask = arrayCenteredMask.GetWidth()*arrayCenteredMask.GetHeight();
	arrayNeighbors.Init(iSizeMask);
	arrayMask.Init(iSizeMask);

	pMask = arrayMask.GetBuffer();
	pCenteredMask = arrayCenteredMask.GetBuffer();
	pNeighbor = arrayNeighbors.GetBuffer();

	for (neighbor.y=-arrayCenteredMask.GetHeight()/2; neighbor.y<=arrayCenteredMask.GetHeight()/2; neighbor.y++)
	{
		for (neighbor.x=-arrayCenteredMask.GetWidth()/2; neighbor.x<=arrayCenteredMask.GetWidth()/2; neighbor.x++)
		{
			*pNeighbor = neighbor;
			*pMask = *pCenteredMask;

			pMask++;
			pCenteredMask++;
			pNeighbor++;
		}
	}
	return Convolve(arrayNeighbors, arrayMask);
}

// Convolve with respect to a given neighborhood
template <class T> template <class S> CArray2D<T> CArray2D<T>::Convolve(const CArray1D<CCouple<int> > &arrayNeighbors, const CArray1D<S> &arrayMask) const
{
	CArray2D<T> arrayRes;
	int iNeighbor;
	const S *pMask;
	int *pOffset;
	const T *pElemY, *pElem;
	T *pElemDestY, *pElemDest;
	CCouple<int> p, p2, minBound, maxBound;
	const CCouple<int> *pNeighbor;
	CArray1D<int> arrayNeighborsOffsets(arrayNeighbors.GetSize());

	// Initialize result and offsets arrays
	arrayRes.Init(iWidth, iHeight);
	arrayNeighborsOffsets.Init(arrayNeighbors.GetSize());

	// Compute bounds
	minBound.Set(numeric_limits<int>::max(), numeric_limits<int>::max());
	maxBound.Set(-numeric_limits<int>::max(), -numeric_limits<int>::max());

	pNeighbor = arrayNeighbors.GetBuffer();
	pOffset = arrayNeighborsOffsets.GetBuffer();
	for (iNeighbor=0;iNeighbor<arrayNeighbors.GetSize();iNeighbor++)
	{
		minBound = coupleMin(minBound, *pNeighbor);
		maxBound = coupleMax(maxBound, *pNeighbor);
		*pOffset = GetOffset(*pNeighbor);

		pNeighbor++;
		pOffset++;
	}

	minBound = -minBound;
	maxBound = CCouple<int>(iWidth-1, iHeight-1)-maxBound;

	// Main loop
	pElemY = this->pElements + minBound.y*iWidth;
	pElemDestY = arrayRes.pElements + minBound.y*iWidth;
	for (p.y=minBound.y; p.y<=maxBound.y; p.y++)
	{
		pElem = pElemY + minBound.x;
		pElemDest = pElemDestY + minBound.x;
		for (p.x=minBound.x; p.x<=maxBound.x; p.x++)
		{
			*pElemDest = 0.0f;
			pOffset = arrayNeighborsOffsets.GetBuffer();
			pMask = arrayMask.GetBuffer();
			for (iNeighbor=0; iNeighbor<arrayNeighbors.GetSize(); iNeighbor++)
			{
				*pElemDest += pElem[*pOffset]*(*pMask);
				pOffset++;
				pMask++;
			}
			pElem++;
			pElemDest++;
		}
		pElemY+=iWidth;
		pElemDestY+=iWidth;
	}

	// Bounds
	CCouple<int> trZero, trSize;
	trZero.Set(0, 0);
	trSize.Set(iWidth-1, iHeight-1);

	// y < minBound.y
	pElem = this->pElements;
	pElemDest = arrayRes.pElements;
	for (p.y=0; p.y<minBound.y; p.y++)
	{
		for (p.x=0; p.x<iWidth; p.x++)
		{
			*pElemDest = 0.0f;
			pNeighbor = arrayNeighbors.GetBuffer();
			pOffset = arrayNeighborsOffsets.GetBuffer();
			pMask = arrayMask.GetBuffer();
			for (iNeighbor=0; iNeighbor<arrayNeighbors.GetSize(); iNeighbor++)
			{
				if ((p+*pNeighbor).IsInRange(trZero, trSize))
					*pElemDest += pElem[*pOffset]*(*pMask);
				else {
					p2 = p+*pNeighbor;
					p2.LimitWithinRange(trZero, trSize);
					*pElemDest += Element(p2)*(*pMask);
				}
				pNeighbor++;
				pOffset++;
				pMask++;
			}
			pElem++;
			pElemDest++;
		}
	}

	// minBound.y <= y <= minBound.y
	pElemY = this->pElements + minBound.y*iWidth;
	pElemDestY = arrayRes.pElements + minBound.y*iWidth;
	for (p.y=minBound.y; p.y<=maxBound.y; p.y++)
	{
		// x < minBound.x
		pElem = pElemY;
		pElemDest = pElemDestY;
		for (p.x=0; p.x<minBound.x; p.x++)
		{
			*pElemDest = 0.0f;
			pNeighbor = arrayNeighbors.GetBuffer();
			pOffset = arrayNeighborsOffsets.GetBuffer();
			pMask = arrayMask.GetBuffer();
			for (iNeighbor=0; iNeighbor<arrayNeighbors.GetSize(); iNeighbor++)
			{
				if ((p+*pNeighbor).IsInRange(trZero, trSize))
					*pElemDest += pElem[*pOffset]*(*pMask);
				else {
					p2 = p+*pNeighbor;
					p2.LimitWithinRange(trZero, trSize);
					*pElemDest += Element(p2)*(*pMask);
				}
				pNeighbor++;
				pOffset++;
				pMask++;
			}
			pElem++;
			pElemDest++;
		}

		// x > maxBound.x
		pElem = pElemY + (maxBound.x+1);
		pElemDest = pElemDestY + (maxBound.x+1);
		for (p.x=maxBound.x+1; p.x<iWidth; p.x++)
		{
			*pElemDest = 0.0f;
			pNeighbor = arrayNeighbors.GetBuffer();
			pOffset = arrayNeighborsOffsets.GetBuffer();
			pMask = arrayMask.GetBuffer();
			for (iNeighbor=0; iNeighbor<arrayNeighbors.GetSize(); iNeighbor++)
			{
				if ((p+*pNeighbor).IsInRange(trZero, trSize))
					*pElemDest += pElem[*pOffset]*(*pMask);
				else {
					p2 = p+*pNeighbor;
					p2.LimitWithinRange(trZero, trSize);
					*pElemDest += Element(p2)*(*pMask);
				}
				pNeighbor++;
				pOffset++;
				pMask++;
			}
			pElem++;
			pElemDest++;
		}

		pElemY+=iWidth;
		pElemDestY+=iWidth;
	}

	// y > maxBound.y
	pElem = this->pElements + (maxBound.y+1)*iWidth;
	pElemDest = arrayRes.pElements + (maxBound.y+1)*iWidth;
	for (p.y=maxBound.y+1; p.y<iHeight; p.y++)
	{
		for (p.x=0; p.x<iWidth; p.x++)
		{
			*pElemDest = 0.0f;
			pNeighbor = arrayNeighbors.GetBuffer();
			pOffset = arrayNeighborsOffsets.GetBuffer();
			pMask = arrayMask.GetBuffer();
			for (iNeighbor=0; iNeighbor<arrayNeighbors.GetSize(); iNeighbor++)
			{
				if ((p+*pNeighbor).IsInRange(trZero, trSize))
					*pElemDest += pElem[*pOffset]*(*pMask);
				else {
					p2 = p+*pNeighbor;
					p2.LimitWithinRange(trZero, trSize);
					*pElemDest += Element(p2)*(*pMask);
				}
				pNeighbor++;
				pOffset++;
				pMask++;
			}
			pElem++;
			pElemDest++;
		}
	}
	return arrayRes;
}

// Convolution with a centered mask and 'marked' array
// An array of booleans of the same size of the current array indicates which elements should be considered in convolution
template <class T> template <class S> CArray2D<T> CArray2D<T>::ConvolveMarked(const CArray2D<S> &arrayCenteredMask, const CArray2D<bool> &arrayMarked) const
{
	CArray1D<CCouple<int> > arrayNeighbors;
	CArray1D<S> arrayMask;
	CCouple<int> neighbor, *pNeighbor;
	S *pMask;
	const S *pCenteredMask;
	int iSizeMask;

	iSizeMask = arrayCenteredMask.GetWidth()*arrayCenteredMask.GetHeight();
	arrayNeighbors.Init(iSizeMask);
	arrayMask.Init(iSizeMask);

	pMask = arrayMask.GetBuffer();
	pCenteredMask = arrayCenteredMask.GetBuffer();
	pNeighbor = arrayNeighbors.GetBuffer();

	for (neighbor.y=-arrayCenteredMask.GetHeight()/2; neighbor.y<=arrayCenteredMask.GetHeight()/2; neighbor.y++)
	{
		for (neighbor.x=-arrayCenteredMask.GetWidth()/2; neighbor.x<=arrayCenteredMask.GetWidth()/2; neighbor.x++)
		{
			*pNeighbor = neighbor;
			*pMask = *pCenteredMask;

			pMask++;
			pCenteredMask++;
			pNeighbor++;
		}
	}
	return ConvolveMarked(arrayNeighbors, arrayMask, arrayMarked);
}

// Convolve with respect to a given neighborhood
// An array of booleans of the same size of the current array indicates which elements should be considered in convolution
template <class T> template <class S> CArray2D<T> CArray2D<T>::ConvolveMarked(const CArray1D<CCouple<int> > &arrayNeighbors, const CArray1D<S> &arrayMask, const CArray2D<bool> &arrayMarked) const
{
	CArray2D<T> arrayRes;
	int iNeighbor;
	const S *pMask;
	S maskSum, maskSumCurrent;

	bool *pMarkedY, *pMarked;
	int *pOffset;
	const T *pElemY, *pElem;
	T *pElemDestY, *pElemDest;
	CCouple<int> p, p2, minBound, maxBound, *pNeighbor;
	CArray1D<int> arrayNeighborsOffsets(arrayNeighbors.GetSize());

    #ifdef ARRAYND_RUNTIME_CHECK
    if (arrayMarked.iWidth!=iWidth || arrayMarked.iHeight!=iHeight)
    {
        cerr<<"ERROR in CArray2D<T>::ConvolveMarked(...): array of booleans has a different size than curren array."<<endl;
        return arrayRes;
    }

    #endif

	// Initialize result and offsets arrays
	arrayRes.Init(iWidth, iHeight);
	arrayNeighborsOffsets.Init(arrayNeighbors.GetSize());

	// Compute bounds
	minBound.Set(numeric_limits<int>::max(), numeric_limits<int>::max());
	maxBound.Set(-numeric_limits<int>::max(), -numeric_limits<int>::max());

	pNeighbor = arrayNeighbors.GetBuffer();
	pOffset = arrayNeighborsOffsets.GetBuffer();
	pMask = arrayMask.GetBuffer();
	maskSum = (S)0;
	for (iNeighbor=0;iNeighbor<arrayNeighbors.GetSize();iNeighbor++)
	{
		minBound = coupleMin(minBound, *pNeighbor);
		maxBound = coupleMax(maxBound, *pNeighbor);
		*pOffset = GetOffset(*pNeighbor);
        maskSum += *pMask;

		pNeighbor++;
		pOffset++;
		pMask++;
	}

	minBound = -minBound;
	maxBound = CCouple<int>(iWidth-1, iHeight-1)-maxBound;

	// Main loop
	pElemY = this->pElements + minBound.y*iWidth;
	pElemDestY = arrayRes.pElements + minBound.y*iWidth;
	pMarkedY = arrayMarked.GetBuffer() + minBound.y*iWidth;

	for (p.y=minBound.y; p.y<=maxBound.y; p.y++)
	{
		pElem = pElemY + minBound.x;
		pElemDest = pElemDestY + minBound.x;
		pMarked = pMarkedY + minBound.x;
		for (p.x=minBound.x; p.x<=maxBound.x; p.x++)
		{
		    if (*pMarked==true)
		    {
                *pElemDest = (T)0;
                maskSumCurrent = (S)0;
                pOffset = arrayNeighborsOffsets.GetBuffer();
                pMask = arrayMask.GetBuffer();
                for (iNeighbor=0; iNeighbor<arrayNeighbors.GetSize(); iNeighbor++)
                {
                    if (pMarked[*pOffset]==true)
                    {
                        *pElemDest += pElem[*pOffset]*(*pMask);
                        maskSumCurrent += *pMask;
                    }

                    pOffset++;
                    pMask++;
                }
                *pElemDest /= maskSumCurrent/maskSum;
		    }
		    else *pElemDest = *pElem;

			pElem++;
			pElemDest++;
			pMarked++;
		}
		pElemY+=iWidth;
		pElemDestY+=iWidth;
		pMarkedY+=iWidth;
	}

	// Bounds
	CCouple<int> trZero, trSize;
	trZero.Set(0, 0);
	trSize.Set(iWidth-1, iHeight-1);

	// y < minBound.y
	pElem = this->pElements;
	pElemDest = arrayRes.pElements;
	for (p.y=0; p.y<minBound.y; p.y++)
	{
		for (p.x=0; p.x<iWidth; p.x++)
		{
			*pElemDest = (T)0;
			pNeighbor = arrayNeighbors.GetBuffer();
			pOffset = arrayNeighborsOffsets.GetBuffer();
			pMask = arrayMask.GetBuffer();
			for (iNeighbor=0; iNeighbor<arrayNeighbors.GetSize(); iNeighbor++)
			{
				if ((p+*pNeighbor).IsInRange(trZero, trSize))
					*pElemDest += pElem[*pOffset]*(*pMask);
				else {
					p2 = p+*pNeighbor;
					p2.LimitWithinRange(trZero, trSize);
					*pElemDest += Element(p2)*(*pMask);
				}
				pNeighbor++;
				pOffset++;
				pMask++;
			}
			pElem++;
			pElemDest++;
		}
	}

	// minBound.y <= y <= minBound.y
	pElemY = this->pElements + minBound.y*iWidth;
	pElemDestY = arrayRes.pElements + minBound.y*iWidth;
	for (p.y=minBound.y; p.y<=maxBound.y; p.y++)
	{
		// x < minBound.x
		pElem = pElemY;
		pElemDest = pElemDestY;
		for (p.x=0; p.x<minBound.x; p.x++)
		{
			*pElemDest = (T)0;
			pNeighbor = arrayNeighbors.GetBuffer();
			pOffset = arrayNeighborsOffsets.GetBuffer();
			pMask = arrayMask.GetBuffer();
			for (iNeighbor=0; iNeighbor<arrayNeighbors.GetSize(); iNeighbor++)
			{
				if ((p+*pNeighbor).IsInRange(trZero, trSize))
					*pElemDest += pElem[*pOffset]*(*pMask);
				else {
					p2 = p+*pNeighbor;
					p2.LimitWithinRange(trZero, trSize);
					*pElemDest += Element(p2)*(*pMask);
				}
				pNeighbor++;
				pOffset++;
				pMask++;
			}
			pElem++;
			pElemDest++;
		}

		// x > maxBound.x
		pElem = pElemY + (maxBound.x+1);
		pElemDest = pElemDestY + (maxBound.x+1);
		for (p.x=maxBound.x+1; p.x<iWidth; p.x++)
		{
			*pElemDest = (T)0;
			pNeighbor = arrayNeighbors.GetBuffer();
			pOffset = arrayNeighborsOffsets.GetBuffer();
			pMask = arrayMask.GetBuffer();
			for (iNeighbor=0; iNeighbor<arrayNeighbors.GetSize(); iNeighbor++)
			{
				if ((p+*pNeighbor).IsInRange(trZero, trSize))
					*pElemDest += pElem[*pOffset]*(*pMask);
				else {
					p2 = p+*pNeighbor;
					p2.LimitWithinRange(trZero, trSize);
					*pElemDest += Element(p2)*(*pMask);
				}
				pNeighbor++;
				pOffset++;
				pMask++;
			}
			pElem++;
			pElemDest++;
		}

		pElemY+=iWidth;
		pElemDestY+=iWidth;
	}

	// y > maxBound.y
	pElem = this->pElements + (maxBound.y+1)*iWidth;
	pElemDest = arrayRes.pElements + (maxBound.y+1)*iWidth;
	for (p.y=maxBound.y+1; p.y<iHeight; p.y++)
	{
		for (p.x=0; p.x<iWidth; p.x++)
		{
			*pElemDest = (T)0;
			pNeighbor = arrayNeighbors.GetBuffer();
			pOffset = arrayNeighborsOffsets.GetBuffer();
			pMask = arrayMask.GetBuffer();
			for (iNeighbor=0; iNeighbor<arrayNeighbors.GetSize(); iNeighbor++)
			{
				if ((p+*pNeighbor).IsInRange(trZero, trSize))
					*pElemDest += pElem[*pOffset]*(*pMask);
				else {
					p2 = p+*pNeighbor;
					p2.LimitWithinRange(trZero, trSize);
					*pElemDest += Element(p2)*(*pMask);
				}
				pNeighbor++;
				pOffset++;
				pMask++;
			}
			pElem++;
			pElemDest++;
		}
	}
	return arrayRes;
}

template <class T> CArray2D<T> CArray2D<T>::SubArray(const CCouple<int> &piStart, const CCouple<int> &piEnd) const
{
    CArray2D<T> arrayDest;
    CCouple<int> piCurrent;
    const T *pElem;
    T *pElemDest;
    int iWidthDest, iHeightDest, iOffsetEndRow;

    #ifdef ARRAYND_RUNTIME_CHECK
    if (piStart.x<0 || piStart.y<0 || piEnd.x>=iWidth || piEnd.y>=iHeight
        || piEnd.x>piStart.x || piEnd.y>piStart.y)
    {
        cerr<<"ERROR in CArray2D<T>::SubArray(...): incorrect bounds. Returning empty array."<<endl;
        return arrayDest;
    }
    #endif

    iWidthDest = piEnd.x-piStart.y+1;
    iHeightDest = piEnd.y-piStart.y+1;
    arrayDest.Init(iWidthDest, iHeightDest);

    iOffsetEndRow = iWidth - iWidthDest;
    pElem = this->pElements + piStart.y*iWidth + piStart.x;
    pElemDest = arrayDest.pElements;

    for (piCurrent.y=piStart.y; piCurrent.y<=piEnd.y; piCurrent.y++)
    {
        for (piCurrent.x=piStart.x; piCurrent.x<=piEnd.x; piCurrent.x++)
        {
            *pElemDest = *pElem;
            pElem++;
            pElemDest++;
        }
        pElem+=iOffsetEndRow;
    }

    return arrayDest;
}

template <class T> class CArray3DIterator
{
  friend class CArray3D<T>;

  // Attributes
  public:
	T *pCurrentElementZ, *pCurrentElementY, *pCurrentElement;
	CTriplet<int> ptCurrent, ptMin, ptMax;
	CArray3D<T> *pArrayParent;

  public:
	inline CArray3DIterator()
	{
		pArrayParent = NULL;
		pCurrentElementZ = pCurrentElementY = pCurrentElement = NULL;
	}
	T *operator ++() // Prefix operator ++
	{
		T *pRet;

		ptCurrent.x++;
		pCurrentElement++;
		if (ptCurrent.x>ptMax.x)
		{
			ptCurrent.x = ptMin.x;
			ptCurrent.y++;

			pCurrentElementY += pArrayParent->iOffsetY;

			if (ptCurrent.y>ptMax.y)
			{
				ptCurrent.y = ptMin.y;
				ptCurrent.z++;

				pCurrentElementZ += pArrayParent->iOffsetZ;
				pCurrentElementY  = pCurrentElementZ + ptMin.y*pArrayParent->iOffsetY;
			}
			pCurrentElement = pCurrentElementY + ptMin.x;
		}
		pRet = pCurrentElement;
		return pRet;
	}
	T *operator ++(int) // Postfix operator ++
	{
		T *pRet;

		pRet = pCurrentElement;

		ptCurrent.x++;
		pCurrentElement++;
		if (ptCurrent.x>ptMax.x)
		{
			ptCurrent.x = ptMin.x;
			ptCurrent.y++;

			pCurrentElementY += pArrayParent->iOffsetY;

			if (ptCurrent.y>ptMax.y)
			{
				ptCurrent.y = ptMin.y;
				ptCurrent.z++;

				pCurrentElementZ += pArrayParent->iOffsetZ;
				pCurrentElementY  = pCurrentElementZ + ptMin.y*pArrayParent->iOffsetY;
			}
			pCurrentElement = pCurrentElementY + ptMin.x;
		}
		return pRet;
	}

	inline CTriplet<int> GetPosition() {return ptCurrent;}
	inline T &Element() {return *pCurrentElement;}
	inline T *ElementPtr() {return pCurrentElement;}
	inline operator T *() {return pCurrentElement;}

	// void Begin() {}
	inline int End() {return (ptCurrent.z>ptMax.z);}
};

template <class T> class CArray3D : protected CArray1D<T>
{
	friend class CArray3DIterator<T>;

  protected:
	int iOffsetY, iOffsetZ;
	int iWidth, iHeight, iDepth;

  // Member functions
  public:
	inline CArray3D();
	inline CArray3D(int, int, int);
	inline CArray3D(const CTriplet<int> &);
	inline CArray3D(const CArray3D<T> &);

	inline ~CArray3D();

	inline int GetWidth() const {return iWidth;}
	inline int GetHeight() const {return iHeight;}
	inline int GetDepth() const {return iDepth;}
	inline CTriplet<int> GetSize() const {return CTriplet<int>(iWidth, iHeight, iDepth);}
	inline int GetOffset(int x, int y, int z) const {return z*iOffsetZ + y*iOffsetY + x;}
	inline int GetOffset(const CTriplet<int> &coord) const {return coord.z*iOffsetZ + coord.y*iOffsetY + coord.x;}

	inline const T *GetBuffer() const {return this->pElements;}
    inline T *GetBuffer() {return this->pElements;}

	bool Init(int, int, int);
	bool Init(const CTriplet<int> &);
	template <class S> bool InitCast(const CArray3D<S> &);
	void Empty();
	void Fill(const T &);

	CArray3D<T> &operator =(const CArray3D<T> &);

	// Access single element (read/write)
	inline T &Element(int x, int y, int z) const
	{
		#ifdef ARRAYND_RUNTIME_CHECK
		if (x<0 || x>=iWidth || y<0 || y>=iHeight || z<0 || z>=iDepth)
		{
			cerr<<"ERROR in CArray3D<T>::Element(int, int, int): accessing element ("<<x<<","<<y<<","<<z<<") out of range [0.."<<iWidth-1<<"]x[0.."<<iHeight-1<<"]x[0.."<<iDepth-1<<"]"<<endl;
			return CArray1D<T>::elementError;
		}
		#endif
		return this->pElements[z*iOffsetZ + y*iOffsetY + x];
	}
	inline T &Element(const CTriplet<int> &coord) const
	{
		#ifdef ARRAYND_RUNTIME_CHECK
		if (coord.x<0 || coord.x>=iWidth || coord.y<0 || coord.y>=iHeight || coord.z<0 || coord.z>=iDepth)
		{
			cerr<<"ERROR in CArray3D<T>::Element(const CTriplet<int> &): accessing element ("<<coord.x<<","<<coord.y<<","<<coord.z<<") out of range [0.."<<iWidth-1<<"]x[0.."<<iHeight-1<<"]x[0.."<<iDepth-1<<"]"<<endl;
			return CArray1D<T>::elementError;
		}
		#endif
		return this->pElements[coord.z*iOffsetZ + coord.y*iOffsetY + coord.x];
	}

	// Linear interpolation between elements (read only)
	inline T GetElementInterpolate(float x, float y, float z) const
	{
		float dx, dy, dz;
		int xi, yi, zi;
		T *pElementTemp;
		T elementInterpolate = (T)0;

		#ifdef ARRAYND_RUNTIME_CHECK
		if (x<0.0f || x>=(float)(iWidth-1) || y<0.0f || y>=(float)(iHeight-1) || z<0.0f || z>=(float)(iDepth-1))
		{
			cerr<<"ERROR in CArray3D<T>::GetElementInterpolate(float,float,float): accessing pixel ("<<x<<","<<y<<") out of range [0.."<<iWidth-1<<"[x[0.."<<iHeight-1<<"[x[0.."<<iDepth-1<<"]"<<endl;
			return elementInterpolate;
		}
		#endif

		xi = (int)x;
		yi = (int)y;
		zi = (int)z;
		dx = x-floor(x);
		dy = y-floor(y);
		dz = z-floor(z);

		// Get address of nearest element with lower integer coordinates
		pElementTemp = this->pElements + zi*this->iOffsetZ + yi*this->iOffsetY + xi;

		elementInterpolate = (T)(
			(1.0f-dx)*(1.0f-dy)*(1.0f-dz) * pElementTemp[0] +
			dx*(1.0f-dy)*(1.0f-dz)        * pElementTemp[1] +
			(1.0f-dx)*dy*(1.0f-dz)        * pElementTemp[iOffsetY] +
			(1.0f-dx)*(1.0f-dy)*dz        * pElementTemp[iOffsetZ] +
			dx*(1.0f-dy)*dz               * pElementTemp[iOffsetZ+1] +
			(1.0f-dx)*dy*dz               * pElementTemp[iOffsetZ+iOffsetY] +
			dx*dy*(1.0f-dz)               * pElementTemp[iOffsetY+1] +
			dx*dy*dz                      * pElementTemp[iOffsetZ+iOffsetY+1]);
		return elementInterpolate;
	}
	inline T GetElementInterpolate(const CTriplet<float> &coord) const {return GetElementInterpolate(coord.x,coord.y,coord.z);}

	// Iterator
	inline CArray3DIterator<T> GetIterator()
	{
		CArray3DIterator<T> iterator;

		iterator.pArrayParent = this;
		iterator.pCurrentElementZ = iterator.pCurrentElementY
			= iterator.pCurrentElement = this->pElements;
		iterator.ptCurrent.init(0,0,0);
		iterator.ptMin.init(0,0,0);
		iterator.ptMax.init(iWidth-1,iHeight-1,iDepth-1);

		return iterator;
	}
	inline CArray3DIterator<T> GetIterator(const CTriplet<int> &ptMin, const CTriplet<int> &ptMax)
	{
		CArray3DIterator<T> iterator;
		CTriplet<int> ptLimitLow(0), ptLimitHigh(GetSize()-CTriplet<int>(1));

		iterator.pArrayParent = this;
		iterator.pCurrentElementZ = this->pElements + ptMin.z*iOffsetZ;
		iterator.pCurrentElementY = iterator.pCurrentElementZ + ptMin.y*iOffsetY;
		iterator.pCurrentElement  = iterator.pCurrentElementY + ptMin.x;
		iterator.ptCurrent = iterator.ptMin = ptMin;
		iterator.ptMax = ptMax;

		if (ptMax.IsInRange(ptLimitLow, ptLimitHigh) && ptMin.IsInRange(ptLimitLow, ptMax))
		{
			iterator.ptCurrent = iterator.ptMin = ptMin;
			iterator.ptMax = ptMax;
		}
		else {
			iterator.ptMax.z = 0;
			iterator.ptCurrent.z = 1;
			cerr<<"ERROR: CArray2D<T>::GetIterator()"<<endl;
		}

		return iterator;
	}

	CArray3D<T> &operator +=(const CArray3D<T> &);
	CArray3D<T> &operator -=(const CArray3D<T> &);
	template <class S> CArray3D<T> &operator *=(const S &);
	template <class S> CArray3D<T> &operator /=(const S &);

	bool operator ==(const CArray3D<T> &) const;
	CArray3D<T> operator +(const CArray3D<T> &) const;
	CArray3D<T> operator -() const;
	CArray3D<T> operator -(const CArray3D<T> &) const;
	template <class S> CArray3D<T> operator *(const S &) const;
	template <class U, class S> friend CArray3D<U> operator *(const S &, const CArray3D<U> &);
	template <class S> CArray3D<T> operator /(const S &) const;

	CArray1D<int> GetNeighborOffsetsConnex6();
	CArray1D<int> GetNeighborOffsetsConnex26();

	T Sum() const {return CArray1D<T>::Sum();}

	// Convolution
	template <class S> CArray3D<T> Convolve(const CArray3D<S> &) const;
	template <class S> CArray3D<T> Convolve(const CArray1D<CCouple<int> > &, const CArray1D<S> &) const;

	// In-place convolution with 1D mask
	template <class S> void ConvolveInPlaceX(const CArray1D<S> &);
	template <class S> void ConvolveInPlaceY(const CArray1D<S> &);
	template <class S> void ConvolveInPlaceZ(const CArray1D<S> &);
};

template <class T> CArray3D<T>::CArray3D():CArray1D<T>()
{
	iWidth = iHeight = iDepth = 0;
	iOffsetY = iOffsetZ = 0;
}

template <class T> CArray3D<T>::CArray3D(int width, int height, int depth):CArray1D<T>(width*height*depth)
{
	iWidth = width;
	iHeight = height;
	iDepth = depth;

	iOffsetY = iWidth;
	iOffsetZ = iWidth*iHeight;
}

template <class T> CArray3D<T>::CArray3D(const CTriplet<int> &iDimensions):CArray1D<T>(iDimensions.x*iDimensions.y*iDimensions.z)
{
	iWidth = iDimensions.x;
	iHeight = iDimensions.y;
	iDepth = iDimensions.z;

	iOffsetY = iWidth;
	iOffsetZ = iWidth*iHeight;
}

template <class T> CArray3D<T>::CArray3D(const CArray3D<T> &array):CArray1D<T>(array)
{
	iWidth = array.iWidth;
	iHeight = array.iHeight;
	iDepth = array.iDepth;

	iOffsetY = iWidth;
	iOffsetZ = iWidth*iHeight;
}

template <class T> CArray3D<T>::~CArray3D()
{
	Empty();
}

template <class T> bool CArray3D<T>::Init(int width, int height, int depth)
{
	iWidth = width;
	iHeight = height;
	iDepth = depth;

	iOffsetY = iWidth;
	iOffsetZ = iWidth*iHeight;

	return CArray1D<T>::Init(iWidth*iHeight*iDepth);
}

template <class T> bool CArray3D<T>::Init(const CTriplet<int> &iDimensions)
{
	iWidth = iDimensions.x;
	iHeight = iDimensions.y;
	iDepth = iDimensions.z;

	iOffsetY = iWidth;
	iOffsetZ = iWidth*iHeight;

	return CArray1D<T>::Init(iWidth*iHeight*iDepth);
}

template <class T> template <class S> bool CArray3D<T>::InitCast(const CArray3D<S> &array)
{
	iWidth = array.GetWidth();
	iHeight = array.GetHeight();
	iDepth = array.GetDepth();

    iOffsetY = iWidth;
	iOffsetZ = iWidth*iHeight;

	if (this->pElements!=NULL)
		delete[] this->pElements;

	this->iSize = iWidth*iHeight*iDepth;

	this->pElements = new T[this->iSize];
	if (this->pElements!=NULL)
	{
		T *pElementsThis = this->pElements;
		S *pElementsParam = array.GetBuffer();
		int i;

		for (i=0;i<this->iSize;i++)
		{
			*pElementsThis = (T)(*pElementsParam);
			pElementsThis++;
			pElementsParam++;
		}
		return true;
	}
	else
		return false;
}

template <class T> void CArray3D<T>::Empty()
{
	iWidth = iHeight = iDepth = 0;
	iOffsetY = iOffsetZ = 0;

	CArray1D<T>::Empty();
}

template <class T> void CArray3D<T>::Fill(const T &t)
{
	CArray1D<T>::Fill(t);
}

template <class T> CArray3D<T> &CArray3D<T>::operator =(const CArray3D<T> &array)
{
	iWidth = array.iWidth;
	iHeight = array.iHeight;
	iDepth = array.iDepth;

	iOffsetY = iWidth;
	iOffsetZ = iWidth*iHeight;

	CArray1D<T>::operator =(array);
	return *this;
}

template <class T> CArray3D<T> &CArray3D<T>::operator +=(const CArray3D<T> &array)
{
	if (iWidth==array.iWidth && iHeight==array.iHeight && iDepth==array.iDepth)
	{
		CArray1D<T>::operator +=(array);
		return *this;
	}
	else
		throw 0;
}

template <class T> CArray3D<T> &CArray3D<T>::operator -=(const CArray3D<T> &array)
{
	if (iWidth==array.iWidth && iHeight==array.iHeight && iDepth==array.iDepth)
	{
		CArray1D<T>::operator -=(array);
		return *this;
	}
	else
		throw 0;
}

template <class T> template <class S> CArray3D<T> &CArray3D<T>::operator *=(const S &t)
{
	CArray1D<T>::operator *=(t);
	return *this;
}

template <class T> template <class S> CArray3D<T> &CArray3D<T>::operator /=(const S &t)
{
	CArray1D<T>::operator /=(t);
	return *this;
}

template <class T> bool CArray3D<T>::operator ==(const CArray3D<T> &array) const
{
	if (iWidth!=array.iWidth || iHeight!=array.iHeight || iDepth!=array.iDepth)
		return false;

	T *pElementsTemp, *pElementsTemp2;
	int i;

	pElementsTemp = this->pElements;
	pElementsTemp2 = array.pElements;
	for (i=0;i<this->iSize && *pElementsTemp==*pElementsTemp2; i++)
	{
		pElementsTemp++;
		pElementsTemp2++;
	}
	if (i==this->iSize) return true;
	else return false;
}

template <class T> CArray3D<T> CArray3D<T>::operator +(const CArray3D<T> &array) const
{
	if (iWidth==array.iWidth && iHeight==array.iHeight && iDepth==array.iDepth)
	{
		CArray3D<T> arrayRes(*this);
		arrayRes.CArray1D<T>::operator +=(array);
		return arrayRes;
	}
	else
		throw 0;
}

template <class T> CArray3D<T> CArray3D<T>::operator -() const
{
	CArray3D<T> arrayResult(iWidth, iHeight, iDepth);
	T *pElementsTemp, *pElementsTempResult;
	int i;

	pElementsTemp = this->pElements;
	pElementsTempResult = arrayResult.pElements;

	for (i=0;i<this->iSize;i++)
	{
		*pElementsTempResult = -(*pElementsTemp);
		pElementsTemp++;
		pElementsTempResult++;
	}
	return arrayResult;
}

template <class T> CArray3D<T> CArray3D<T>::operator -(const CArray3D<T> &array) const
{
	if (iWidth==array.iWidth && iHeight==array.iHeight && iDepth==array.iDepth)
	{
		CArray3D<T> arrayRes(*this);
		arrayRes.CArray1D<T>::operator -=(array);
		return arrayRes;
	}
	else
		throw 0;
}

template <class T> template <class S> CArray3D<T> CArray3D<T>::operator *(const S &t) const
{
	CArray3D<T> arrayRes(*this);
	arrayRes.CArray1D<T>::operator *=(t);
	return arrayRes;
}

template <class T, class S> CArray3D<T> operator *(const S &t, const CArray3D<T> &array)
{
	CArray3D<T> arrayRes(array);
	arrayRes.CArray1D<T>::operator *=(t);
	return arrayRes;
}

template <class T> template <class S> CArray3D<T> CArray3D<T>::operator /(const S &t) const
{
	CArray3D<T> arrayRes(*this);
	arrayRes.CArray1D<T>::operator /=(t);
	return arrayRes;
}

template <class T>  CArray1D<int> CArray3D<T>::GetNeighborOffsetsConnex6()
{
	CArray1D<int> arrayOffsets(6);

	arrayOffsets[0] = -1;
	arrayOffsets[1] = 1;
	arrayOffsets[2] = -iOffsetY;
	arrayOffsets[3] = iOffsetY;
	arrayOffsets[4] = -iOffsetZ;
	arrayOffsets[5] = iOffsetZ;

	return arrayOffsets;
}

template <class T>  CArray1D<int> CArray3D<T>::GetNeighborOffsetsConnex26()
{
	CArray1D<int> arrayOffsets(26);
	CTriplet<int> pi;
	int i;

	i = 0;
	for (pi.z=-1;pi.z<=1;pi.z++)
	{
		for (pi.y=-1;pi.y<=1;pi.y++)
		{
			for (pi.x=-1;pi.x<=1;pi.x++)
			{
				if (pi.x!=0 || pi.y!=0 || pi.z!=0)
				{
					arrayOffsets[i++] = pi.z*iOffsetZ + pi.y*iWidth + pi.x;
				}
			}
		}
	}
	return arrayOffsets;
}

// Convolution with a centered mask
template <class T> template <class S> CArray3D<T> CArray3D<T>::Convolve(const CArray3D<S> &arrayCenteredMask) const
{
	CArray1D<CTriplet<int> > arrayNeighbors;
	CArray1D<S> arrayMask;
	CTriplet<int> neighbor, *pNeighbor;
	S *pMask, *pCenteredMask;
	int iSizeMask;

	iSizeMask = arrayCenteredMask.GetWidth()*arrayCenteredMask.GetHeight()*arrayCenteredMask.GetDepth();
	arrayNeighbors.Init(iSizeMask);
	arrayMask.Init(iSizeMask);

	pMask = arrayMask.GetBuffer();
	pCenteredMask = arrayCenteredMask.GetBuffer();
	pNeighbor = arrayNeighbors.GetBuffer();
	for (neighbor.z=-arrayCenteredMask.GetDepth()/2; neighbor.z<=arrayCenteredMask.GetDepth()/2; neighbor.z++)
	{
		for (neighbor.y=-arrayCenteredMask.GetHeight()/2; neighbor.y<=arrayCenteredMask.GetHeight()/2; neighbor.y++)
		{
			for (neighbor.x=-arrayCenteredMask.GetWidth()/2; neighbor.x<=arrayCenteredMask.GetWidth()/2; neighbor.x++)
			{
				*pNeighbor = neighbor;
				*pMask = *pCenteredMask;

				pMask++;
				pCenteredMask++;
				pNeighbor++;
			}
		}
	}

	return Convolve(arrayNeighbors, arrayMask);
}

// Convolve with respect to a given neighborhood
template <class T> template <class S> CArray3D<T> CArray3D<T>::Convolve(const CArray1D<CCouple<int> > &arrayNeighbors, const CArray1D<S> &arrayMask) const
{
	CArray3D<T> arrayRes;
	int iNeighbor;
	S *pMask;
	int *pOffset;
	T *pElemZ, *pElemY, *pElem;
	T *pElemDestZ, *pElemDestY, *pElemDest;
	CTriplet<int> p, p2, minBound, maxBound, *pNeighbor;
	CArray1D<int> arrayNeighborsOffsets(arrayNeighbors.GetSize());

	// Initialize result and offsets arrays
	arrayRes.Init(iWidth, iHeight, iDepth);
	arrayNeighborsOffsets.Init(arrayNeighbors.GetSize());

	// Compute bounds
	minBound.Set(numeric_limits<int>::max(), numeric_limits<int>::max(), numeric_limits<int>::max());
	maxBound.Set(-numeric_limits<int>::max(), -numeric_limits<int>::max(), -numeric_limits<int>::max());

	pNeighbor = arrayNeighbors.GetBuffer();
	pOffset = arrayNeighborsOffsets.GetBuffer();
	for (iNeighbor=0;iNeighbor<arrayNeighbors.GetSize();iNeighbor++)
	{
		minBound = tripletMin(minBound, *pNeighbor);
		maxBound = tripletMax(maxBound, *pNeighbor);
		*pOffset = GetOffset(*pNeighbor);

		pNeighbor++;
		pOffset++;
	}

	minBound = -minBound;
	maxBound = CTriplet<int>(iWidth-1, iHeight-1, iDepth-1)-maxBound;

	// Main loop
	pElemZ = this->pElements + minBound.z*iOffsetZ;
	pElemDestZ = arrayRes.pElements + minBound.z*iOffsetZ;
	for (p.z=minBound.z; p.z<=maxBound.z; p.z++)
	{
		pElemY = pElemZ + minBound.y*iOffsetY;
		pElemDestY = pElemDestZ + minBound.y*iOffsetY;
		for (p.y=minBound.y; p.y<=maxBound.y; p.y++)
		{
			pElem = pElemY + minBound.x;
			pElemDest = pElemDestY + minBound.x;
			for (p.x=minBound.x; p.x<=maxBound.x; p.x++)
			{
				*pElemDest = 0.0f;
				pOffset = arrayNeighborsOffsets.GetBuffer();
				pMask = arrayMask.GetBuffer();
				for (iNeighbor=0; iNeighbor<arrayNeighbors.GetSize(); iNeighbor++)
				{
					*pElemDest += pElem[*pOffset]*(*pMask);
					pOffset++;
					pMask++;
				}
				pElem++;
				pElemDest++;
			}
			pElemY+=iOffsetY;
			pElemDestY+=iOffsetY;
		}
		pElemZ+=iOffsetZ;
		pElemDestZ+=iOffsetZ;
	}

	// Bounds
	CTriplet<int> trZero, trSize;
	trZero.Set(0, 0, 0);
	trSize.Set(iWidth-1, iHeight-1, iDepth-1);

	// z < minBound.z
	pElem = this->pElements;
	pElemDest = arrayRes.pElements;
	for (p.z=0; p.z<minBound.z; p.z++)
	{
		for (p.y=0; p.y<iHeight; p.y++)
		{
			for (p.x=0; p.x<iWidth; p.x++)
			{
				*pElemDest = 0.0f;
				pNeighbor = arrayNeighbors.GetBuffer();
				pOffset = arrayNeighborsOffsets.GetBuffer();
				pMask = arrayMask.GetBuffer();
				for (iNeighbor=0; iNeighbor<arrayNeighbors.GetSize(); iNeighbor++)
				{
					if ((p+*pNeighbor).IsInRange(trZero, trSize))
						*pElemDest += pElem[*pOffset]*(*pMask);
					else {
						p2 = p+*pNeighbor;
						p2.LimitWithinRange(trZero, trSize);
						*pElemDest += Element(p2)*(*pMask);
					}
					pNeighbor++;
					pOffset++;
					pMask++;
				}
				pElem++;
				pElemDest++;
			}
		}
	}

	// minBound.z <= z <=maxBound.z
	pElemZ = this->pElements + minBound.z*iOffsetZ;
	pElemDestZ = arrayRes.pElements + minBound.z*iOffsetZ;
	for (p.z=minBound.z; p.z<=maxBound.z; p.z++)
	{
		// y < minBound.y
		pElem = pElemZ;
		pElemDest = pElemDestZ;
		for (p.y=0; p.y<minBound.y; p.y++)
		{
			for (p.x=0; p.x<iWidth; p.x++)
			{
				*pElemDest = 0.0f;
				pNeighbor = arrayNeighbors.GetBuffer();
				pOffset = arrayNeighborsOffsets.GetBuffer();
				pMask = arrayMask.GetBuffer();
				for (iNeighbor=0; iNeighbor<arrayNeighbors.GetSize(); iNeighbor++)
				{
					if ((p+*pNeighbor).IsInRange(trZero, trSize))
						*pElemDest += pElem[*pOffset]*(*pMask);
					else {
						p2 = p+*pNeighbor;
						p2.LimitWithinRange(trZero, trSize);
						*pElemDest += Element(p2)*(*pMask);
					}
					pNeighbor++;
					pOffset++;
					pMask++;
				}
				pElem++;
				pElemDest++;
			}
		}

		// minBound.y <= y <= minBound.y
		pElemY = pElemZ + minBound.y*iOffsetY;
		pElemDestY = pElemDestZ + minBound.y*iOffsetY;
		for (p.y=minBound.y; p.y<=maxBound.y; p.y++)
		{
			// x < minBound.x
			pElem = pElemY;
			pElemDest = pElemDestY;
			for (p.x=0; p.x<minBound.x; p.x++)
			{
				*pElemDest = 0.0f;
				pNeighbor = arrayNeighbors.GetBuffer();
				pOffset = arrayNeighborsOffsets.GetBuffer();
				pMask = arrayMask.GetBuffer();
				for (iNeighbor=0; iNeighbor<arrayNeighbors.GetSize(); iNeighbor++)
				{
					if ((p+*pNeighbor).IsInRange(trZero, trSize))
						*pElemDest += pElem[*pOffset]*(*pMask);
					else {
						p2 = p+*pNeighbor;
						p2.LimitWithinRange(trZero, trSize);
						*pElemDest += Element(p2)*(*pMask);
					}
					pNeighbor++;
					pOffset++;
					pMask++;
				}
				pElem++;
				pElemDest++;
			}

			// x > maxBound.x
			pElem = pElemY + (maxBound.x+1);
			pElemDest = pElemDestY + (maxBound.x+1);
			for (p.x=maxBound.x+1; p.x<iWidth; p.x++)
			{
				*pElemDest = 0.0f;
				pNeighbor = arrayNeighbors.GetBuffer();
				pOffset = arrayNeighborsOffsets.GetBuffer();
				pMask = arrayMask.GetBuffer();
				for (iNeighbor=0; iNeighbor<arrayNeighbors.GetSize(); iNeighbor++)
				{
					if ((p+*pNeighbor).IsInRange(trZero, trSize))
						*pElemDest += pElem[*pOffset]*(*pMask);
					else {
						p2 = p+*pNeighbor;
						p2.LimitWithinRange(trZero, trSize);
						*pElemDest += Element(p2)*(*pMask);
					}
					pNeighbor++;
					pOffset++;
					pMask++;
				}
				pElem++;
				pElemDest++;
			}

			pElemY+=iOffsetY;
			pElemDestY+=iOffsetY;
		}

		// y > maxBound.y
		pElem = pElemZ + (maxBound.y+1)*iOffsetY;
		pElemDest = pElemDestZ + (maxBound.y+1)*iOffsetY;
		for (p.y=maxBound.y+1; p.y<iHeight; p.y++)
		{
			for (p.x=0; p.x<iWidth; p.x++)
			{
				*pElemDest = 0.0f;
				pNeighbor = arrayNeighbors.GetBuffer();
				pOffset = arrayNeighborsOffsets.GetBuffer();
				pMask = arrayMask.GetBuffer();
				for (iNeighbor=0; iNeighbor<arrayNeighbors.GetSize(); iNeighbor++)
				{
					if ((p+*pNeighbor).IsInRange(trZero, trSize))
						*pElemDest += pElem[*pOffset]*(*pMask);
					else {
						p2 = p+*pNeighbor;
						p2.LimitWithinRange(trZero, trSize);
						*pElemDest += Element(p2)*(*pMask);
					}
					pNeighbor++;
					pOffset++;
					pMask++;
				}
				pElem++;
				pElemDest++;
			}
		}

		pElemZ+=iOffsetZ;
		pElemDestZ+=iOffsetZ;
	}

	// z > maxBound.z
	pElem = this->pElements + (maxBound.z+1)*iOffsetZ;
	pElemDest = arrayRes.pElements + (maxBound.z+1)*iOffsetZ;
	for (p.z=maxBound.z+1; p.z<iDepth; p.z++)
	{
		for (p.y=0; p.y<iHeight; p.y++)
		{
			for (p.x=0; p.x<iWidth; p.x++)
			{
				*pElemDest = 0.0f;
				pNeighbor = arrayNeighbors.GetBuffer();
				pOffset = arrayNeighborsOffsets.GetBuffer();
				pMask = arrayMask.GetBuffer();
				for (iNeighbor=0; iNeighbor<arrayNeighbors.GetSize(); iNeighbor++)
				{
					if ((p+*pNeighbor).IsInRange(trZero, trSize))
						*pElemDest += pElem[*pOffset]*(*pMask);
					else {
						p2 = p+*pNeighbor;
						p2.LimitWithinRange(trZero, trSize);
						*pElemDest += Element(p2)*(*pMask);
					}
					pNeighbor++;
					pOffset++;
					pMask++;
				}
				pElem++;
				pElemDest++;
			}
		}
	}
	return arrayRes;
}


// In-place convolution with 1D mask in x-dimension
template <class T> template <class S> void CArray3D<T>::ConvolveInPlaceX(const CArray1D<S> &arrayCenteredRowMask)
{
	CArray1D<T> arrayOutputRowTemp;
	CTriplet<int> p;
	int iNeighbor, iMaskHalfSize;
	T *pElem, *pRowStart, *pNeighbor;
	T *pOutputRowTemp;
	const S *pMask;
	int iStartX, iEndX;

	// Initialize temporary row array
	arrayOutputRowTemp.Init(iWidth);

	iMaskHalfSize = arrayCenteredRowMask.GetSize()/2;
	iStartX = iMaskHalfSize;
	iEndX = iWidth - 1 - iMaskHalfSize;

	pRowStart = this->pElements + iStartX;
	for (p.z=0; p.z<iDepth; p.z++)
	{
		for (p.y=0; p.y<iHeight; p.y++)
		{
			pElem = pRowStart;
			pOutputRowTemp = arrayOutputRowTemp.GetBuffer() + iStartX;
			for (p.x=iStartX; p.x<=iEndX; p.x++)
			{
				*pOutputRowTemp = (T)0;
				pNeighbor = pElem - iMaskHalfSize;
				pMask = arrayCenteredRowMask.GetBuffer();
				for (iNeighbor=0; iNeighbor<arrayCenteredRowMask.GetSize(); iNeighbor++)
				{
					*pOutputRowTemp += (*pNeighbor)*(*pMask);
					pNeighbor++;
					pMask++;
				}
				pElem++;
				pOutputRowTemp++;
			}

			pElem = pRowStart;
			pOutputRowTemp = arrayOutputRowTemp.GetBuffer() + iStartX;
			for (p.x=iStartX; p.x<=iEndX; p.x++)
			{
				*pElem = *pOutputRowTemp;
				pElem++;
				pOutputRowTemp++;
			}

			pRowStart+=iWidth;
		}
	}
}

// In-place convolution with 1D mask in y-dimension
template <class T> template <class S> void CArray3D<T>::ConvolveInPlaceY(const CArray1D<S> &arrayCenteredColumnMask)
{
	CArray1D<T> arrayOutputColumnTemp;
	CTriplet<int> p;
	int iNeighbor, iMaskHalfSize;
	T *pElem, *pColumnStart, *pNeighbor;
	T *pOutputColumnTemp;
	const S *pMask;
	int iStartY, iEndY;

	// Initialize temporary row array
	arrayOutputColumnTemp.Init(iHeight);

	iMaskHalfSize = arrayCenteredColumnMask.GetSize()/2;
	iStartY = iMaskHalfSize;
	iEndY = iHeight - 1 - iMaskHalfSize;

	pColumnStart = this->pElements + iStartY*iWidth;
	for (p.z=0; p.z<iDepth; p.z++)
	{
		for (p.x=0; p.x<iWidth; p.x++)
		{
			pElem = pColumnStart;
			pOutputColumnTemp = arrayOutputColumnTemp.GetBuffer() + iStartY;

			for (p.y=iStartY; p.y<=iEndY; p.y++)
			{
				*pOutputColumnTemp = (T)0;
				pNeighbor = pElem - iMaskHalfSize*iWidth;
				pMask = arrayCenteredColumnMask.GetBuffer();
				for (iNeighbor=0; iNeighbor<arrayCenteredColumnMask.GetSize(); iNeighbor++)
				{
					*pOutputColumnTemp += (*pNeighbor)*(*pMask);
					pNeighbor+=iWidth;
					pMask++;
				}
				pElem+=iWidth;
				pOutputColumnTemp++;
			}

			pElem = pColumnStart;
			pOutputColumnTemp = arrayOutputColumnTemp.GetBuffer() + iStartY;
			for (p.y=iStartY; p.y<=iEndY; p.y++)
			{
				*pElem = *pOutputColumnTemp;
				pElem+=iWidth;
				pOutputColumnTemp++;
			}

			pColumnStart++;
		}
		pColumnStart += iOffsetZ - iWidth;
	}
}

// In-place convolution with 1D mask in z-dimension
template <class T> template <class S> void CArray3D<T>::ConvolveInPlaceZ(const CArray1D<S> &arrayCenteredRowZMask)
{
	CArray1D<T> arrayOutputRowZTemp;
	CTriplet<int> p;
	int iNeighbor, iMaskHalfSize;
	T *pElem, *pRowZStart, *pNeighbor;
	T *pOutputRowZTemp;
	const S *pMask;
	int iStartZ, iEndZ;

	// Initialize temporary row array
	arrayOutputRowZTemp.Init(iDepth);

	iMaskHalfSize = arrayCenteredRowZMask.GetSize()/2;
	iStartZ = iMaskHalfSize;
	iEndZ = iDepth - 1 - iMaskHalfSize;

	pRowZStart = this->pElements + iStartZ*iOffsetZ;
	for (p.y=0; p.y<iHeight; p.y++)
	{
		for (p.x=0; p.x<iWidth; p.x++)
		{
			pElem = pRowZStart;
			pOutputRowZTemp = arrayOutputRowZTemp.GetBuffer() + iStartZ;
			for (p.z=iStartZ; p.z<=iEndZ; p.z++)
			{
				*pOutputRowZTemp = (T)0;
				pNeighbor = pElem - iMaskHalfSize*iOffsetZ;
				pMask = arrayCenteredRowZMask.GetBuffer();
				for (iNeighbor=0; iNeighbor<arrayCenteredRowZMask.GetSize(); iNeighbor++)
				{
					*pOutputRowZTemp += (*pNeighbor)*(*pMask);
					pNeighbor+=iOffsetZ;
					pMask++;
				}
				pElem+=iOffsetZ;
				pOutputRowZTemp++;
			}

			pElem = pRowZStart;
			pOutputRowZTemp = arrayOutputRowZTemp.GetBuffer() + iStartZ;
			for (p.z=iStartZ; p.z<=iEndZ; p.z++)
			{
				*pElem = *pOutputRowZTemp;
				pElem+=iOffsetZ;
				pOutputRowZTemp++;
			}

			pRowZStart++;
		}
	}
}

#endif
