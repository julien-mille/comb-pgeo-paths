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
path.h

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

#ifndef _PATH_H_
#define _PATH_H_

#include <list>

#include "couple.h"
#include "image2d.h"
#include "arraynd.h"

using namespace std;

class CIntegerPath2D;
class CRealPath2D;

class CIntegerPath2D : public list<CCouple<int> >
{
  friend class CRealPath2D;

  // Member variables
  public:
	CImage2DByteRGBPixel rgbVertices, rgbEdges;
	bool bDisplayVertices, bDisplayEdges;
    float fSortingValue;

  // Member functions
  public:
	CIntegerPath2D()
	{
		rgbVertices.Set(0,255,0);
		rgbEdges.Set(0,0,0);
		bDisplayVertices = true;
		bDisplayEdges = true;
	}
	CIntegerPath2D(const CRealPath2D &);

    // CIntegerPath2D is used to implement the discretized 4-connected admissible paths
    // Before searching combination of admissible paths, we sort them in admissible
    // sets with respect to some value (exteriority). We need to overload comparison
    // operators to allow the STL sort() funtion
    bool operator <(const CIntegerPath2D &path) const {
        return fSortingValue<path.fSortingValue;
    }
    bool operator >(const CIntegerPath2D &path) const {
        return fSortingValue>path.fSortingValue;
    }

	float GetLength() const;

	// Compute area with discrete Green's theorem (should be a digital 4-connected closed curve)
	float GetArea4connected() const;

    void SetDigitalLineSegment4connected(const CCouple<int> &, const CCouple<int> &);
    void AddDigitalPath4connected(const CIntegerPath2D &);

	void DrawInImage(CImage2D &, int) const;

	// Draw the path by considering points as pointels (pixel corners) instead of pixels
	void DrawInImagePixelCorners(CImage2D &, int) const;
};

class CRealPath2D : public list<CCouple<float> >
{
  friend class CIntegerPath2D;

  // Member variables
  public:
	CImage2DByteRGBPixel rgbVertices, rgbEdges;
	bool bDisplayVertices, bDisplayEdges;
    float fSortingValue;

  // Member functions
  public:
	CRealPath2D()
	{
		rgbVertices.Set(0,255,0);
		rgbEdges.Set(0,0,0);
		bDisplayVertices = true;
		bDisplayEdges = true;
	}
	CRealPath2D(const CIntegerPath2D &);

    // Admissible paths are implemented in curve with real point coordinates
    // Before searching combination of admissible paths, we sort them in admissible
    // sets with respect to some value (exteriority). We need to overload comparison
    // operators to allow the STL sort() funtion
    bool operator <(const CRealPath2D &path) const {
        return fSortingValue<path.fSortingValue;
    }
    bool operator >(const CRealPath2D &path) const {
        return fSortingValue>path.fSortingValue;
    }

	float GetLength() const;
	void InitLine(const CCouple<float> &, const CCouple<float> &);
	void LaplacianSmooth();

	void Translate(const CCouple<float> &);

    float GetSignedAreaWithLineSegment() const;

    CRealPath2D Resampled(float) const;

	void DrawInImage(CImage2D &, int) const;
};

#endif
