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

#ifndef _ARRAYNDFLOAT_H_
#define _ARRAYNDFLOAT_H_

#include "arraynd.h"
#include "triplet.h"

// Possible schemes for 1st-order finite differences
#define ARRAYNDFLOAT_BACKWARD 0
#define ARRAYNDFLOAT_FORWARD  1
#define ARRAYNDFLOAT_CENTERED 2 // The centered method is applied by default

template <> class CArray1D<float>
{
  // Static members
#ifdef ARRAYND_RUNTIME_CHECK
  protected:
	// Element returned in case of invalid access with operator [] or function GetElementInterpolate()
	static float elementError;
#endif

  // Members
  protected:
	int iSize;
	float *pElements;

  // Methods
  public:
	inline CArray1D() {iSize = 0; pElements = NULL;}
	CArray1D(int);
	CArray1D(const CArray1D<float> &);

	~CArray1D();

	inline int GetSize() const {return iSize;}
	inline const float *GetBuffer() const {return pElements;}
    inline float *GetBuffer() {return pElements;}

	bool Init(int);
	bool Init(int, float *);

	// MSVC linker returns unresolved external error if not defined inline
	template <class S> bool InitCast(const CArray1D<S> &array)
	{
		if (pElements!=NULL)
			delete[] pElements;

		iSize = array.GetSize();

		pElements = new float[iSize];
		if (pElements!=NULL)
		{
			float *pElementsThis = pElements;
			const S *pElementsParam = array.GetBuffer();
			int i;

			for (i=0;i<iSize;i++)
			{
				*pElementsThis = (float)(*pElementsParam);
				pElementsThis++;
				pElementsParam++;
			}
			return true;
		}
		else
			return false;
	}

	void Empty();
	void Fill(const float &);

	CArray1D<float> &operator =(const CArray1D<float> &);

	// Access single element (read/write)
	inline float &operator [](int index) const
	{
		#ifdef ARRAYND_RUNTIME_CHECK
		if (index<0 || index>=iSize)
		{
			cerr<<"ERROR in CArray1D<float>::operator [](int): accessing element at index "<<index<<" out of range [0.."<<iSize-1<<"]"<<endl;
			return elementError;
		}
		#endif
		return pElements[index];
	}

	// Linear interpolation between elements (read only)
	inline float GetElementInterpolate(float fIndex) const
	{
		float fDiff;
		int iIndex;
		float elementInterpolate;

		#ifdef ARRAYND_RUNTIME_CHECK
		if (fIndex<0.0f || fIndex>=(float)(iSize-1))
		{
			cerr<<"ERROR in CArray1D<float>::GetElementInterpolate(float): accessing element at index "<<fIndex<<" out of range [0.."<<iSize-1<<"["<<endl;
			return elementError;
		}
		#endif

		iIndex = (int)fIndex;
		fDiff = fIndex-floor(fIndex);
		elementInterpolate = (1.0f-fDiff)*pElements[iIndex] + fDiff*pElements[iIndex+1];

		return elementInterpolate;
	}

	CArray1D<float> &operator +=(const CArray1D<float> &);
	CArray1D<float> &operator -=(const CArray1D<float> &);
	CArray1D<float> &operator *=(const float &);
	CArray1D<float> &operator /=(const float &);

	CArray1D<float> operator +(const CArray1D<float> &) const;
	CArray1D<float> operator -() const;
	CArray1D<float> operator -(const CArray1D<float> &) const;
	CArray1D<float> operator *(const float &) const;
	friend CArray1D<float> operator *(const float &, const CArray1D<float> &);
	friend CArray1D<float> MultiplyElements(const CArray1D<float> &, const CArray1D<float> &);

	CArray1D<float> operator /(const float &) const;

	CArray1D<float> &AbsElements();
	CArray1D<float> &SquareElements();
	CArray1D<float> &SqrtElements();

	// Norms and distances
	CCouple<float> GetMinMax() const;
	float Sum() const;
	float L1Norm() const;
	float L2Norm() const;
	float InfiniteNorm() const;

	friend float L1Distance(const CArray1D<float> &, const CArray1D<float> &);
	friend float L2Distance(const CArray1D<float> &, const CArray1D<float> &);
	friend float InfiniteDistance(const CArray1D<float> &, const CArray1D<float> &);

	void Limit(float, float);
	void Normalize(float, float, float &, float &);

    // Convolution
	CArray1D<float> Convolve(const CArray1D<float> &) const;
	CArray1D<float> Convolve(const CArray1D<int> &, const CArray1D<float> &) const;

    // Faster convolution with centered mask (borders are ignored)
    void ConvolveNoBorder(const CArray1D<float> &);

	// Gaussian and related kernel
	void SetGaussianKernel(float, int halfSize=0);
};

template <> class CArray2D<float> : protected CArray1D<float>
{
  protected:
	int iWidth, iHeight;

  // Methods
  public:
	inline CArray2D():CArray1D<float>() {iWidth = iHeight = 0;}
	CArray2D(int, int);
	CArray2D(const CCouple<int> &);
	CArray2D(const CArray2D<float> &);

	~CArray2D();

	inline int GetWidth() const {return iWidth;}
	inline int GetHeight() const {return iHeight;}
	inline CCouple<int> GetSize() const {return CCouple<int>(iWidth, iHeight);}
	inline int GetOffset(int x, int y) const {return y*iWidth+x;}
	inline int GetOffset(const CCouple<int> &coord) const {return coord.y*iWidth+coord.x;}

	inline const float *GetBuffer() const {return pElements;}
    inline float *GetBuffer() {return pElements;}

	bool Init(int, int);
	bool Init(int, int, float *);
	bool Init(const CCouple<int> &);

	// MSVC linker returns unresolved external error if not defined inline
	template <class S> bool InitCast(const CArray2D<S> &array)
	{
		if (pElements!=NULL)
			delete[] pElements;

		iWidth = array.GetWidth();
		iHeight = array.GetHeight();

		iSize = iWidth*iHeight;

		pElements = new float[iSize];
		if (pElements!=NULL)
		{
			float *pElementsThis = pElements;
			const S *pElementsParam = array.GetBuffer();
			int i;

			for (i=0;i<iSize;i++)
			{
				*pElementsThis = (float)(*pElementsParam);
				pElementsThis++;
				pElementsParam++;
			}
			return true;
		}
		else
			return false;
	}

	void Empty();
	void Fill(float);
	void FillBoundary(float);

	CArray2D<float> &operator =(const CArray2D<float> &);

	// Access single element (read/write)
	// Access single element (read/write)
	inline float &Element(int x, int y) const
	{
		#ifdef ARRAYND_RUNTIME_CHECK
		if (x<0 || x>=iWidth || y<0 || y>=iHeight)
		{
			cerr<<"ERROR in CArray2D<float>::Element(int, int): accessing element ("<<x<<","<<y<<") out of range [0.."<<iWidth-1<<"]x[0.."<<iHeight-1<<"]"<<endl;
			return CArray1D<float>::elementError;
		}
		#endif
		return pElements[y*iWidth+x];
	}
	inline float &Element(const CCouple<int> &coord) const
	{
		#ifdef ARRAYND_RUNTIME_CHECK
		if (coord.x<0 || coord.x>=iWidth || coord.y<0 || coord.y>=iHeight)
		{
			cerr<<"ERROR in CArray2D<float>::Element(const CCouple<int> &): accessing element ("<<coord.x<<","<<coord.y<<") out of range [0.."<<iWidth-1<<"]x[0.."<<iHeight-1<<"]"<<endl;
			return CArray1D<float>::elementError;
		}
		#endif
		return pElements[coord.y*iWidth+coord.x];
	}

	// Linear interpolation between elements (read only)
	inline float GetElementInterpolate(float x, float y) const
	{
		float dx, dy;
		int xi, yi;
		const float *pElementTemp;
		float elementInterpolate;

		#ifdef ARRAYND_RUNTIME_CHECK
		if (x<0.0f || x>=(float)(iWidth-1) || y<0.0f || y>=(float)(iHeight-1))
		{
			cerr<<"ERROR in CArray2D<float>::GetElementInterpolate(float,float): accessing pixel ("<<x<<","<<y<<") out of range [0.."<<iWidth-1<<"[x[0.."<<iHeight-1<<"["<<endl;
			return 0.0f;
		}
		#endif

		xi = (int)x;
		yi = (int)y;
		dx = x-floor(x);
		dy = y-floor(y);

		// Get address of nearest element with lower integer coordinates
		pElementTemp = this->pElements + yi*this->iWidth + xi;

		elementInterpolate =
			(1.0f-dx)*(1.0f-dy) * pElementTemp[0] +
			dx*(1.0f-dy)        * pElementTemp[1] +
			(1.0f-dx)*dy        * pElementTemp[iWidth] +
			dx*dy               * pElementTemp[iWidth+1];
		return elementInterpolate;
	}
	inline float GetElementInterpolate(const CCouple<float> &coord) const {return GetElementInterpolate(coord.x,coord.y);}

	CArray2D<float> &operator +=(const CArray2D<float> &);
	CArray2D<float> &operator -=(const CArray2D<float> &);
	CArray2D<float> &operator *=(float);
	CArray2D<float> &operator /=(float);

	CArray2D<float> operator +(const CArray2D<float> &) const;
	CArray2D<float> operator -() const;
	CArray2D<float> operator -(const CArray2D<float> &) const;
	CArray2D<float> operator *(float) const;
	friend CArray2D<float> operator *(float, const CArray2D<float> &);
	friend CArray2D<float> MultiplyElements(const CArray2D<float> &, const CArray2D<float> &);
	CArray2D<float> operator /(float) const;

	CArray2D<float> &AbsElements();
	CArray2D<float> &SquareElements();
	CArray2D<float> &SqrtElements();

	// Norms and distances
	CCouple<float> GetMinMax() const;
	float Sum() const;
	float L1Norm() const;
	float L2Norm() const;
	float InfiniteNorm() const;

	friend float L1Distance(const CArray2D<float> &, const CArray2D<float> &);
	friend float L2Distance(const CArray2D<float> &, const CArray2D<float> &);
	friend float InfiniteDistance(const CArray2D<float> &, const CArray2D<float> &);

	CArray1D<int> GetNeighborOffsetsConnex4() const;
	CArray1D<int> GetNeighborOffsetsConnex8() const;

	void Limit(float, float);
	void Normalize(float, float, float &, float &);

	// Differentiation
	CArray2D<float> DerivativeX(int iOrder=1, int iScheme=ARRAYNDFLOAT_CENTERED) const;
	CArray2D<float> DerivativeY(int iOrder=1, int iScheme=ARRAYNDFLOAT_CENTERED) const;
	CArray2D<CCouple<float> > Gradient(int iScheme=ARRAYNDFLOAT_CENTERED) const;
	CArray2D<float> GradientNorm() const;
	CArray2D<float> Laplacian() const;

	void SetL2DistanceMap(const CCouple<int> &);

	// Convolution
	CArray2D<float> Convolve(const CArray2D<float> &) const;
	CArray2D<float> Convolve(const CArray1D<CCouple<int> > &, const CArray1D<float> &) const;

	// Gaussian and related kernels
	void SetGaussianKernel(float, int halfSize=0);
	void SetLaplacianOfGaussianKernel(float, int halfSize=0);

	CArray2D<float> SubArray(const CCouple<int> &, const CCouple<int> &) const;
};

#endif
