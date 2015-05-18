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
voronoigraph.h

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

#ifndef _VORONOIGRAPH_H_
#define _VORONOIGRAPH_H_

#include <list>
#include <vector>
#include "arrayndfloat.h"
#include "path.h"

using namespace std;

// Voronoi edge points are pointels (pixel corners) which have neighboring pixels
// belonging to different Voronoi cells.
// A Voronoi edge point may be a saddle point, i.e. a local minima of the combined action map
// along a Voronoi edge (medial curve), in which case the normal vector to the medial curve
// will be estimated
class CVoronoiEdgePoint
{
  // Member variables
  public:
    CCouple<int> piPos; // Pointel (pixel corner)
    float fAction; // Interpolated action
    bool bLocalMinimum;
    CCouple<float> vfNormal; // Normal to the medial curve (if it is a saddle point)

  // Member functions
  public:
    CVoronoiEdgePoint() {}

    // Saddle points are sorted with respect to their action value so we need to overload
    // comparison operators to allow the STL sort() function
    bool operator <(const CVoronoiEdgePoint &saddle) const {return fAction<saddle.fAction;}
    bool operator >(const CVoronoiEdgePoint &saddle) const {return fAction>saddle.fAction;}
    bool operator <=(const CVoronoiEdgePoint &saddle) const {return fAction<=saddle.fAction;}
    bool operator >=(const CVoronoiEdgePoint &saddle) const {return fAction>=saddle.fAction;}
    bool operator ==(const CVoronoiEdgePoint &saddle) const {return fAction==saddle.fAction;}
};

// Voronoi edges (medial curves) are represented by digital 4-connected curves, made up of
// Voronoi edge points, located on pointels. Two successive Voronoi edge points define a segment
// between two adjacent pixels, which belong to two different Voronoi cells
// See Fig. 3 in appendix B.1
class CVoronoiEdge : public list<CVoronoiEdgePoint>
{
  public:
    void push_back(const CCouple<int> &);
	void push_front(const CCouple<int> &);

    // Laplacian smoothing of action values
    void SmoothActionsAlongEdge();

    // Compute saddle points as second-order local minima of action values
    // along a Voronoi edge. Normal vectors to the medial curve are estimated
    // Saddle points are stored in the output list
    // Params: list of saddle points [out]
	void ComputeLocalMinimaAlongEdge(list<CVoronoiEdgePoint> &) const;

    CIntegerPath2D ToIntegerPath() const;
};

// The Voronoi graph is made up of the action (geodesic distance) map and the array of Voronoi cells,
// which contains integer labels corresponding to source points of the action map
class CVoronoiGraph
{
  // Member variables
  protected:
	// Action and Voronoi maps
	// These arrays should be initialized by a geodesic-distance propagation algorithm
	// before any Voronoi edge or set of saddle points can be extracted
	const CArray2D<float> *pArrayAction;
    const CArray2D<int> *pArrayVoronoi;

  // Member functions
  public:
    CVoronoiGraph()
    {
        pArrayAction = NULL;
        pArrayVoronoi = NULL;
    }

    // Initialize action and Voronoi maps
    // Arrays should have the same size
    // Params: pointer to action map [in], pointer to Voronoi map [in]
    // Return value : TRUE if arrays have the same size, and at least 5x5, FALSE otherwise
    bool SetActionAndVoronoiMaps(const CArray2D<float> *, const CArray2D<int> *);

    // Find the Voronoi edge point having the lowest interpolated action value
    // Params: the point located on a Voronoi edge [out] with the smallest action, if any
    // Return value: TRUE if a Voronoi edge point was found, FALSE otherwise
    bool GetVoronoiEdgePointMinimalAction(CCouple<int> &) const;

    // Compute action at pointel using bilinear interpolation of action values
    // of the 4 neighboring pixels
    // Params: any point within the Voronoi map, except border [in]
    // Return value : interpolated action
    float GetActionAtPointel(const CCouple<int> &) const;

    // Build Voronoi edge by edge-linking procedure
    // The starting point should be a valid Voronoi edge point
    // The obtained digital 4-connected curve may be closed or open, regarding the
    // topology of the Voronoi map
    // Params: point used as a starting position for the edge-linking algorithm [in],
    //         Voronoi edge (curve with integer coordinates) [out]
	void GetEdgeFromPointel(const CCouple<int> &, CVoronoiEdge &) const;
};

#endif
