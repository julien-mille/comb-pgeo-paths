/*
Copyright 2015 Julien Mille

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

/*
histogram.h

Header file of library which implements the algorithm combining
piecewise-geodesic paths to build closed contours, described in

[MBC15] J. Mille, S. Bougleux and L. Cohen. Combination of piecewise-geodesic
        paths for interactive segmentation. International Journal of Computer
        Vision, 112(1):1-22, 2015.

All equation, section and algorithm numbers in the comments refer to this paper
or in the supplemental document

Preprint versions of the paper and the supplemental document may be found at
http://liris.cnrs.fr/~jmille/doc/ijcv15temp.pdf
http://liris.cnrs.fr/~jmille/doc/ijcv15supp.pdf

If you use this code for research purposes, please cite the aforementioned paper
in any resulting publication.
*/

#ifndef _HISTOGRAM_H_
#define _HISTOGRAM_H_

#include "arrayndfloat.h"
#include <vector>

class CHistogram1D : protected CArray1D<float>
{
  // Member variables
  protected:
    // Number of bins
    int iNbBins;

    // Gaussian kernel for smoothing
	float fGaussianStdDeviation;
	CArray1D<float> arrayGaussian;

	int iPadding;

    // Intensity range of unscaled input values
	float fRangeMin, fRangeMax;

	float fSumWeights;

  // Member functions
  public:
	CHistogram1D();

	virtual void SetNbBins(int);
	virtual void SetGaussianStdDeviation(float);
    virtual void SetRange(float, float);

    virtual float GetSumWeights() const {return fSumWeights;}

	virtual void Empty()
	{
        Fill(0.0f);
        fSumWeights = 0.0f;
    }

	virtual void Destroy()
	{
        iNbBins = 0;
        iPadding = 0;
        CArray1D<float>::Empty();
    }

	// Add and subtract intensity values (in the original unscaled grayscale space)
	virtual void AddElement(float, float fWeight=1.0f);
	virtual void SubtractElement(float, float fWeight=1.0f);

	// Add and subtract scaled intensity values
	// If scaled inetnsity values are precomputed, these are much faster
	// than AddElement(...) and SubtractElement(...)
	virtual void AddElementScaled(int, float fWeight=1.0f);
	virtual void SubtractElementScaled(int, float fWeight=1.0f);

    // Add and subtract histograms
    virtual CHistogram1D &operator +=(const CHistogram1D &);
    virtual CHistogram1D &operator -=(const CHistogram1D &);
    virtual CHistogram1D operator +(const CHistogram1D &) const;
    virtual CHistogram1D operator -(const CHistogram1D &) const;

	// Get bin value (with linear interpolation)
	virtual float GetValue(float) const;

    // Create an image of integer scaled intensities, where
    // values are scaled and rounded according to the range of the histogram
    // Filling histograms will be much faster if handling scaled values
    virtual void MakeImageScaled(const CArray2D<float> &, CArray2D<int> &) const;

    // In-place convolution with Gaussian filter
	virtual void GaussianSmooth();

     // Compute BhattacharyyaCoefficient between two probability distributions functions
    // (PDFs are obtained by normalizing bin values by the sum of weights)
    virtual float BhattacharyyaCoefficient(const CHistogram1D &) const;
};

class CHistogram3D : protected CArray3D<float>
{
  // Member variables
  protected:
    // Number of bins per component
    CTriplet<int> iNbBins;

    // Gaussian kernel for smoothing
	float fGaussianStdDeviation;
	CArray1D<float> arrayGaussian;

	int iPadding;

    // 3D color range of unscaled input values
	CTriplet<float> fRangeMin, fRangeMax;

	float fSumWeights;

  // Member functions
  public:
	CHistogram3D();

	virtual void SetNbBins(const CTriplet<int> &);
	virtual void SetGaussianStdDeviation(float);
    virtual void SetRange(const CTriplet<float> &, const CTriplet<float> &);

    virtual float GetSumWeights() const {return fSumWeights;}

	virtual void Empty()
	{
        Fill(0.0f);
        fSumWeights = 0.0f;
    }

	virtual void Destroy()
	{
        iNbBins.Set(0,0,0);
        iPadding = 0;
        CArray3D<float>::Empty();
    }

	// Add and subtract colors (in the original unscaled color space)
	virtual void AddElement(const CTriplet<float> &, float fWeight=1.0f);
	virtual void SubtractElement(const CTriplet<float> &, float fWeight=1.0f);

	// Add and subtract scaled colors
	// If scaled colors are precomputed, these are much faster
	// than AddElement(...) and SubtractElement(...)
	virtual void AddElementScaled(const CTriplet<int> &, float fWeight=1.0f);
	virtual void SubtractElementScaled(const CTriplet<int> &, float fWeight=1.0f);

    // Add and subtract histograms
    virtual CHistogram3D &operator +=(const CHistogram3D &);
    virtual CHistogram3D &operator -=(const CHistogram3D &);
    virtual CHistogram3D operator +(const CHistogram3D &) const;
    virtual CHistogram3D operator -(const CHistogram3D &) const;

	// Get bin value (with linear interpolation)
	virtual float GetValue(const CTriplet<float> &) const;

    // Create an image of integer scaled colors, where
    // values are scaled and rounded according to the range of the histogram
    // Filling histograms will be much faster if handling scaled values
    virtual void MakeImageScaled(const CArray2D<CTriplet<float> > &, CArray2D<CTriplet<int> > &) const;

    // In-place convolution with separable Gaussian filter
	virtual void GaussianSmooth();

    // Compute BhattacharyyaCoefficient between two probability distributions functions
    // (PDFs are obtained by normalizing bin values by the sum of weights)
    virtual float BhattacharyyaCoefficient(const CHistogram3D &) const;
};

#endif
