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
voronoigraph.cpp

Source file of library which implements the algorithm combining
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

#include <limits>

#include "voronoigraph.h"

void CVoronoiEdge::push_back(const CCouple<int> &piPos)
{
    CVoronoiEdgePoint vorPtNew;

    vorPtNew.piPos = piPos;
    vorPtNew.fAction = 0.0f;
    vorPtNew.bLocalMinimum = false;

    list<CVoronoiEdgePoint>::push_back(vorPtNew);
}

void CVoronoiEdge::push_front(const CCouple<int> &piPos)
{
    CVoronoiEdgePoint vorPtNew;

    vorPtNew.piPos = piPos;
    vorPtNew.fAction = 0.0f;
    vorPtNew.bLocalMinimum = false;

    list<CVoronoiEdgePoint>::push_front(vorPtNew);
}

void CVoronoiEdge::SmoothActionsAlongEdge()
{
	CArray1D<float> vectNewActions;
	CVoronoiEdge::iterator itPoint, itPointPrev, itPointNext;
	unsigned int iPoint;

	if (size()<3)
		return;

	vectNewActions.Init(size()-2);

	for (itPoint=++begin(), iPoint=0; itPoint!=--end(); itPoint++, iPoint++)
	{
	    itPointPrev = itPoint;
	    itPointPrev--;
	    itPointNext = itPoint;
	    itPointNext++;

		vectNewActions[iPoint] = (itPoint->fAction + (itPointPrev->fAction + itPointNext->fAction)*0.5f)*0.5f;
	}

	for (itPoint=++begin(), iPoint=0; itPoint!=--end(); itPoint++, iPoint++)
		itPoint->fAction = vectNewActions[iPoint];
}


void CVoronoiEdge::ComputeLocalMinimaAlongEdge(list<CVoronoiEdgePoint> &listSaddles) const
{
    CVoronoiEdge::const_iterator itPoint, itPointPrev, itPointPrev2, itPointNext, itPointNext2;
    CVoronoiEdgePoint vorPtSaddle;

    listSaddles.clear();

    for (itPoint=++(++begin()); itPoint!=--(--end()); itPoint++)
    {
        itPointPrev = itPoint;
        itPointPrev--;
        itPointPrev2 = itPointPrev;
        itPointPrev2--;

        itPointNext = itPoint;
        itPointNext++;
        itPointNext2 = itPointNext;
        itPointNext2++;

        if (itPointPrev2->fAction > itPointPrev->fAction &&
            itPointPrev->fAction > itPoint->fAction &&
            itPoint->fAction < itPointNext->fAction &&
            itPointNext->fAction < itPointNext2->fAction)
        {
            vorPtSaddle = *itPoint;
            vorPtSaddle.vfNormal = ((CCouple<float>)((itPointNext->piPos - itPointPrev->piPos).Perp())).Normalized();
            listSaddles.push_back(vorPtSaddle);
        }
    }

    listSaddles.sort();
}

CIntegerPath2D CVoronoiEdge::ToIntegerPath() const
{
    CVoronoiEdge::const_iterator itVoronoiPoint;
    CIntegerPath2D pathVoronoi;

    for (itVoronoiPoint=begin(); itVoronoiPoint!=end(); itVoronoiPoint++)
        pathVoronoi.push_back(itVoronoiPoint->piPos);

    return pathVoronoi;
}

bool CVoronoiGraph::SetActionAndVoronoiMaps(const CArray2D<float> *pAct, const CArray2D<int> *pVor)
{
    if (pAct->GetSize()!=pVor->GetSize())
    {
        cerr<<"ERROR in CVoronoiGraph::SetActionAndVoronoiMaps(...): arrays have different sizes"<<endl;
        return false;
    }

    if (pAct->GetWidth()<5 || pAct->GetHeight()<5)
    {
        cerr<<"ERROR in CVoronoiGraph::SetActionAndVoronoiMaps(...): arrays are too small.";
        cerr<<" They should be at least 5x5 so that Voronoi edges can be extracted"<<endl;
        return false;
    }

    pArrayAction = pAct;
    pArrayVoronoi = pVor;

    return true;
}

bool CVoronoiGraph::GetVoronoiEdgePointMinimalAction(CCouple<int> &piPointelMinDist) const
{
    CArray1D<bool> arraySurroundingCells;
    int i, iSize, iNeighbor, iVoronoiMin, iVoronoiMax, iNbSurroundingCells;
    CCouple<int> piSize, piCurrent;
    const int *pVoronoi, *pVoronoiNeighbor;
    float fActionMin;
    bool bFoundEdgePoint; // Indicate whether at least one Voronoi edge point was found, with a valid action

    piSize = pArrayVoronoi->GetSize();

    iSize = piSize.x*piSize.y;
    pVoronoi = pArrayVoronoi->GetBuffer();
    iVoronoiMin = *pVoronoi;
    iVoronoiMax = *pVoronoi;
    for (i=0; i<iSize; i++)
    {
        if (*pVoronoi<iVoronoiMin)
            iVoronoiMin = *pVoronoi;
        else if (*pVoronoi>iVoronoiMax)
            iVoronoiMax = *pVoronoi;
        pVoronoi++;
    }

    arraySurroundingCells.Init(iVoronoiMax-iVoronoiMin+1);

    int arrayNeighborsOffsets[4] = {-piSize.x-1, -piSize.x, -1, 0};

    piPointelMinDist.Set(-1,-1);
    fActionMin = numeric_limits<float>::max();

    bFoundEdgePoint = false;
    pVoronoi = pArrayVoronoi->GetBuffer() + pArrayVoronoi->GetOffset(2,2);
    for (piCurrent.y=2; piCurrent.y<piSize.y-1; piCurrent.y++)
    {
        for (piCurrent.x=2; piCurrent.x<piSize.x-1; piCurrent.x++)
        {
            arraySurroundingCells.Fill(false);
            iNbSurroundingCells = 0;
            for (iNeighbor=0; iNeighbor<4; iNeighbor++)
            {
                pVoronoiNeighbor = pVoronoi + arrayNeighborsOffsets[iNeighbor];
                if (arraySurroundingCells[*pVoronoiNeighbor-iVoronoiMin]==false)
                {
                    iNbSurroundingCells++;
                    arraySurroundingCells[*pVoronoiNeighbor-iVoronoiMin] = true;
                }
            }

            if (iNbSurroundingCells>=2)
            {
                if (pArrayAction->Element(piCurrent)<fActionMin)
                {
                    fActionMin = pArrayAction->Element(piCurrent);
                    piPointelMinDist = piCurrent;
                    bFoundEdgePoint = true;
                }
            }
            pVoronoi++;
        }
        pVoronoi+=3;
    }

    return bFoundEdgePoint;
}

float CVoronoiGraph::GetActionAtPointel(const CCouple<int> &piPointel) const
{
    float fAction;
    const float *pAction;
    int iWidth = pArrayAction->GetWidth();

    pAction = pArrayAction->GetBuffer() + pArrayAction->GetOffset(piPointel);
    fAction = (*pAction + pAction[-1] + pAction[-iWidth] + pAction[-iWidth-1])*0.25f;

    return fAction;
}

void CVoronoiGraph::GetEdgeFromPointel(const CCouple<int> &pointelStart, CVoronoiEdge &pathVoronoi) const
{
	CCouple<int> piNeighbor, piMin, piMax;
	CCouple<CCouple<int> > pairPixel;
	CCouple<int> pointelCurrent, pointelNeighbor;
	CArray2D<bool> arrayVisited;
	bool bEnd = false, bNeighborFound;
	unsigned int iNeighbor;

    vector<CCouple<int> > vectNeighbors;
    vector<CCouple<CCouple<int> > > vectNeighborPixelPairs;

    // Up
    vectNeighbors.push_back(CCouple<int>(0,-1));
    pairPixel.x.Set(-1,-1);
    pairPixel.y.Set(0,-1);
    vectNeighborPixelPairs.push_back(pairPixel);

    // Right
    vectNeighbors.push_back(CCouple<int>(1,0));
    pairPixel.x.Set(0,-1);
    pairPixel.y.Set(0,0);
    vectNeighborPixelPairs.push_back(pairPixel);

    // Down
    vectNeighbors.push_back(CCouple<int>(0,1));
    pairPixel.x.Set(-1,0);
    pairPixel.y.Set(0,0);
    vectNeighborPixelPairs.push_back(pairPixel);

    // Left
    vectNeighbors.push_back(CCouple<int>(-1,0));
    pairPixel.x.Set(-1,-1);
    pairPixel.y.Set(-1,0);
    vectNeighborPixelPairs.push_back(pairPixel);

    piMin.Set(3,3);
    piMax = pArrayVoronoi->GetSize()-CCouple<int>(4,4);

    arrayVisited.Init(pArrayVoronoi->GetSize()+CCouple<int>(1,1));
    arrayVisited.Fill(false);

	pathVoronoi.clear();

	pointelCurrent = pointelStart;
	arrayVisited.Element(pointelCurrent) = true;
	pathVoronoi.push_back(pointelCurrent);
    pathVoronoi.back().fAction = GetActionAtPointel(pointelCurrent);

	while (bEnd==false)
	{
	    iNeighbor = 0;
	    bNeighborFound = false;
	    while (iNeighbor<vectNeighbors.size() && bNeighborFound==false)
	    {
	        pointelNeighbor = pointelCurrent + vectNeighbors[iNeighbor];

	        if (arrayVisited.Element(pointelNeighbor)==false &&
                pArrayVoronoi->Element(pointelCurrent + vectNeighborPixelPairs[iNeighbor].x)!=
                pArrayVoronoi->Element(pointelCurrent + vectNeighborPixelPairs[iNeighbor].y))
                bNeighborFound = true;
            else
                iNeighbor++;
	    }

	    if (bNeighborFound==true)
        {
            pointelCurrent += vectNeighbors[iNeighbor];
            arrayVisited.Element(pointelCurrent) = true;
            pathVoronoi.push_back(pointelCurrent);
            pathVoronoi.back().fAction = GetActionAtPointel(pointelCurrent);

            if (pointelCurrent.IsInRange(piMin, piMax)==false)
                bEnd = true;
        }
        else bEnd = true;
	}

    pointelCurrent = pointelStart;
    bEnd = false;
	while (bEnd==false)
	{
	    iNeighbor = 0;
	    bNeighborFound = false;
	    while (iNeighbor<vectNeighbors.size() && bNeighborFound==false)
	    {
	        pointelNeighbor = pointelCurrent + vectNeighbors[iNeighbor];

	        if (arrayVisited.Element(pointelNeighbor)==false &&
                pArrayVoronoi->Element(pointelCurrent + vectNeighborPixelPairs[iNeighbor].x)!=
                pArrayVoronoi->Element(pointelCurrent + vectNeighborPixelPairs[iNeighbor].y))
                bNeighborFound = true;
            else
                iNeighbor++;
	    }

	    if (bNeighborFound==true)
        {
            pointelCurrent += vectNeighbors[iNeighbor];
            arrayVisited.Element(pointelCurrent) = true;
            pathVoronoi.push_front(pointelCurrent);
            pathVoronoi.front().fAction = GetActionAtPointel(pointelCurrent);
            if (pointelCurrent.IsInRange(piMin, piMax)==false)
                bEnd = true;
        }
        else bEnd = true;
	}
}
