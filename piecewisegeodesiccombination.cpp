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
piecewisegeodesiccombination.cpp

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

#include <iostream>

#include "geocomp/include/fm2.hh"
#include "piecewisegeodesiccombination.h"
#include "voronoigraph.h"
// #include "graphicwnd.h"


/**********************************************************************
*                CPiecewiseGeodesicCombinationBase                  *
**********************************************************************/

CPiecewiseGeodesicCombinationBase::CPiecewiseGeodesicCombinationBase()
{
    // Set the algorithm to find the best combination of paths
    iSearchType = SEARCHTYPE_GREEDY;

    // Weights used in potential
	fWeightSmoothnessPotential = 0.01f;
	fWeightGradientPotential = 15.0f;
	fGradientScale = 1.0f;

	// Weights of energies
	fWeightSimplicity = 2.0f;
	fWeightEdge = 1.0f;
	fWeightRegion = 2.0f;

    iNbAdmissiblePathsMax = 5;

    iNbBinsPerComponent = 32;

	// Display parameters
	bDisplayVertices = true;
	bDisplayAllAdmissiblePaths = false;
	bDisplayVoronoiEdges = false;

	rgbEdgeColor.Set(255, 0, 0);
	rgbVertexColor.Set(255, 255, 0);
	iVertexWidth = 3;

	pInputImage = NULL;
}

void CPiecewiseGeodesicCombinationBase::AddVertex(const CCouple<float> &pfNew)
{
    CVertex2D vertexNew;

    listListsAdmissiblePaths.clear();
    listListsAdmissiblePathsDigital.clear();
    listVoronoiEdges.clear();

    vertexNew.pfPos = pfNew;
    listVertices.push_back(vertexNew);
}

void CPiecewiseGeodesicCombinationBase::DrawInImageRGB(CImage2D &im, int iZoom) const
{
	list<list<CRealPath2D> >::const_iterator itListPaths;
    list<CRealPath2D>::const_iterator itPath;
    unsigned int idxListPaths, idxPath;

	if (im.GetBitsPerPixel()!=24 && im.GetBitsPerPixel()!=32)
	{
		cerr<<"ERROR in CPiecewiseGeodesicCombinationBase::DrawInImageRGB(...): image is not RGB"<<endl;
		return;
	}

    // Draw admissible paths
    if (listVertices.size()==2 && listListsAdmissiblePaths.size()<=2)
    {
        itListPaths = listListsAdmissiblePaths.begin();
        for (itPath=itListPaths->begin(), idxPath=0; itPath!=itListPaths->end(); itPath++, idxPath++)
            if (bDisplayAllAdmissiblePaths==true || idxPath==vectBestCombination[0])
                itPath->DrawInImage(im, iZoom);
    }
    else {
        for (itListPaths=listListsAdmissiblePaths.begin(), idxListPaths=0;
            itListPaths!=listListsAdmissiblePaths.end(); itListPaths++, idxListPaths++)
            for (itPath=itListPaths->begin(), idxPath=0; itPath!=itListPaths->end(); itPath++, idxPath++)
                if (bDisplayAllAdmissiblePaths==true || idxPath==vectBestCombination[idxListPaths])
                    itPath->DrawInImage(im, iZoom);
    }

    // Draw Voronoi edges
    if (bDisplayVoronoiEdges==true)
    {
        list<CVoronoiEdge>::const_iterator itVoronoiEdge;
        CIntegerPath2D pathVoronoi;

        for (itVoronoiEdge=listVoronoiEdges.begin();itVoronoiEdge!=listVoronoiEdges.end();itVoronoiEdge++)
        {
            pathVoronoi = itVoronoiEdge->ToIntegerPath();
            pathVoronoi.rgbEdges.Set(0,0,0);
            pathVoronoi.DrawInImage(im, iZoom);
        }
    }

	// Draw vertices
	if (bDisplayVertices==true)
	{
        list<CVertex2D>::const_iterator itVertex;

        // Compute integer coordinates of vertices with current zooming factor
        for (itVertex=listVertices.begin(); itVertex!=listVertices.end(); itVertex++)
            im.DrawFilledCircleRGB(
                (CCouple<int>)(itVertex->pfPos*(float)iZoom) + CCouple<int>(iZoom/2, iZoom/2),
                iVertexWidth, rgbVertexColor);;
	}
}

void CPiecewiseGeodesicCombinationBase::MakeBinaryMask(CImage2D &imgMask) const
{
	list<list<CIntegerPath2D> >::const_iterator itListPathsDigital;
    list<CIntegerPath2D>::const_iterator itPathDigital;
    CIntegerPath2D::const_iterator itPoint;
    unsigned int idxListPaths, idxPath;

	if (imgMask.GetSize()!=pInputImage->GetSize() || imgMask.GetBitsPerPixel()!=8)
        imgMask.Create(pInputImage->GetWidth(), pInputImage->GetHeight(), 8);

    imgMask.Clear(255);

    for (itListPathsDigital=listListsAdmissiblePathsDigital.begin(), idxListPaths=0;
        itListPathsDigital!=listListsAdmissiblePathsDigital.end(); itListPathsDigital++, idxListPaths++)
    {
        itPathDigital = itListPathsDigital->begin();
        idxPath = 0;
        while (idxPath!=vectBestCombination[idxListPaths])
        {
            itPathDigital++;
            idxPath++;
        }

        for (itPoint=itPathDigital->begin(); itPoint!=itPathDigital->end(); itPoint++)
            imgMask.Pixel(*itPoint) = 254;
    }

    imgMask.FloodFill(CCouple<int>(0,0), 0);
}

void CPiecewiseGeodesicCombinationBase::Empty()
{
	listVertices.clear();
    listListsAdmissiblePaths.clear();
    listListsAdmissiblePathsDigital.clear();

    listVoronoiEdges.clear();
}

bool CPiecewiseGeodesicCombinationBase::UpdatePotential()
{
	int i, iSize;
	float *pGradient, *pPotential;

	if (arrayPotential.GetSize()!=arrayGradientNorm.GetSize())
        arrayPotential.Init(arrayGradientNorm.GetSize());

	iSize = arrayPotential.GetWidth()*arrayPotential.GetHeight();
    pPotential = arrayPotential.GetBuffer();
    pGradient = arrayGradientNorm.GetBuffer();

    for (i=0; i<iSize; i++)
    {
        *pPotential = fWeightSmoothnessPotential + max(0.0f, 1.0f-fWeightGradientPotential*(*pGradient));
        pPotential++;
        pGradient++;
    }

	return true;
}

void CPiecewiseGeodesicCombinationBase::UpdateAllPaths()
{
	list<CVertex2D>::const_iterator itVertex, itVertexNext, itVertexEnd;
    list<list<CRealPath2D> >::iterator itListPaths;
    list<CRealPath2D>::iterator itPath;
	unsigned int idxVertex;

	listListsAdmissiblePaths.clear();
	listListsAdmissiblePathsDigital.clear();
	listVoronoiEdges.clear();

    if (listVertices.size()<2)
    {
        cerr<<"ERROR in CPiecewiseGeodesicCombinationBase::UpdateAllPaths(): Need at least 2 vertices"<<endl;
        return;
    }

    if (iNbAdmissiblePathsMax<1)
    {
        cerr<<"ERROR in CPiecewiseGeodesicCombinationBase::UpdateAllPaths(): Need at least 1 admissible path per pair of successive vertices"<<endl;
        return;
    }

    if (listVertices.size()==2 && iNbAdmissiblePathsMax==1)
        cout<<"WARNING in CPiecewiseGeodesicCombinationBase::UpdateAllPaths(): Only 2 vertices and the shortest path between them. Cannot build closed contour"<<endl;

    for (itVertex=listVertices.begin(), idxVertex=0; idxVertex<listVertices.size()-1; itVertex++, idxVertex++)
    {
        itVertexNext = itVertex;
        itVertexNext++;

        listListsAdmissiblePaths.push_back(list<CRealPath2D>());

        if (iNbAdmissiblePathsMax==1)
        {
            listListsAdmissiblePaths.back().push_back(CRealPath2D());
            ComputeShortestPath(itVertex->pfPos, itVertexNext->pfPos, listListsAdmissiblePaths.back().back());
            listListsAdmissiblePaths.back().back().rgbEdges = rgbEdgeColor;
        }
        else {
            listVoronoiEdges.push_back(CVoronoiEdge());
            ComputeAdmissiblePaths(itVertex->pfPos, itVertexNext->pfPos, listListsAdmissiblePaths.back(), listVoronoiEdges.back());
        }
    }

    if (listVertices.size()>2)
    {
        itVertex = --listVertices.end();
        itVertexNext = listVertices.begin();

        listListsAdmissiblePaths.push_back(list<CRealPath2D>());

        if (iNbAdmissiblePathsMax==1)
        {
            listListsAdmissiblePaths.back().push_back(CRealPath2D());
            ComputeShortestPath(itVertex->pfPos, itVertexNext->pfPos, listListsAdmissiblePaths.back().back());
            listListsAdmissiblePaths.back().back().rgbEdges = rgbEdgeColor;
        }
        else {
            listVoronoiEdges.push_back(CVoronoiEdge());
            ComputeAdmissiblePaths(itVertex->pfPos, itVertexNext->pfPos, listListsAdmissiblePaths.back(), listVoronoiEdges.back());
        }
    }
    else {
        // Copy paths from the first admissible set into the second one
        // and reverse them
        list<CRealPath2D> *pListFirst, *pListSecond;

        listListsAdmissiblePaths.push_back(list<CRealPath2D>());

        pListFirst = &(listListsAdmissiblePaths.front());
        pListSecond = &(listListsAdmissiblePaths.back());

        for (itPath=pListFirst->begin(); itPath!=pListFirst->end(); itPath++)
        {
            pListSecond->push_back(*itPath);
            pListSecond->back().reverse();
        }

        // Duplicate Voronoi edge
        if (iNbAdmissiblePathsMax>1)
            listVoronoiEdges.push_back(listVoronoiEdges.front());
    }


    // Convert to digital paths
    list<CIntegerPath2D> *pListPathsDigital;

    for (itListPaths=listListsAdmissiblePaths.begin(); itListPaths!=listListsAdmissiblePaths.end(); itListPaths++)
    {
        listListsAdmissiblePathsDigital.push_back(list<CIntegerPath2D>());
        pListPathsDigital = &(listListsAdmissiblePathsDigital.back());
        for (itPath=itListPaths->begin(); itPath!=itListPaths->end(); itPath++)
        {
            itPath->bDisplayVertices = false;
            itPath->bDisplayEdges = true;

            pListPathsDigital->push_back(CIntegerPath2D());
            ConvertRealPathToDigitalPathSimple(*itPath, pListPathsDigital->back());
        }
    }

    if (iNbAdmissiblePathsMax>1)
    {
        SortAdmissiblePathsByExteriority();

        if (listVertices.size()>2)
        {
            if (iSearchType==SEARCHTYPE_GREEDY)
                FindBestCombinationGreedy();
            else if (iSearchType==SEARCHTYPE_EXHAUSTIVE)
                FindBestCombinationExhaustive();
        }
        else
            FindBestCombinationTwoVertices();
    }
    else {
        vectBestCombination.resize(listVertices.size());
        for (idxVertex=0; idxVertex<vectBestCombination.size(); idxVertex++)
            vectBestCombination[idxVertex] = 0;
    }

    // Set display properties of paths
    for (idxVertex=0, itListPaths=listListsAdmissiblePaths.begin();
        idxVertex<listListsAdmissiblePaths.size(); idxVertex++, itListPaths++)
    {
        for (itPath=itListPaths->begin(); itPath!=itListPaths->end(); itPath++)
        {
            itPath->bDisplayVertices = false;
            itPath->bDisplayEdges = true;
        }
    }
}

void CPiecewiseGeodesicCombinationBase::ComputeShortestPath(const CCouple<float> &pfStart, const CCouple<float> &pfEnd, CRealPath2D &realPath)
{
	int iWidth, iHeight;

	iWidth = pInputImage->GetWidth();
	iHeight = pInputImage->GetHeight();

    GeoComp::Grid8<int> grid(iHeight, iWidth);
    GeoComp::FM2<float, GeoComp::Grid8<int> > fastMarching(grid);

    // Action map and its first-order derivatives
    CArray2D<float> arrayAction, arrayActionDerX, arrayActionDerY;
    float *U, *DUx, *DUy, *P;

    arrayAction.Init(iWidth, iHeight);
    arrayActionDerX.Init(iWidth, iHeight);
    arrayActionDerY.Init(iWidth, iHeight);

    U = arrayAction.GetBuffer();
    DUx = arrayActionDerX.GetBuffer();
    DUy = arrayActionDerY.GetBuffer();

    P = arrayPotential.GetBuffer();

    // Propagate action map using Fast Marching
    fastMarching.distance(grid.index((int)(pfStart.y+0.5f), (int)(pfStart.x+0.5f)),P,U,DUy,DUx);

    GradientDescentActionMap(pfEnd, arrayActionDerX, arrayActionDerY, realPath);
    realPath.reverse();
}

void CPiecewiseGeodesicCombinationBase::ComputeAdmissiblePaths(const CCouple<float> &pfStart, const CCouple<float> &pfEnd, list<CRealPath2D> &listPaths, CVoronoiEdge &pathVoronoi)
{
    CImage2D im2d;
	int iWidth, iHeight;
	CVoronoiGraph voronoiGraph;

	iWidth = pInputImage->GetWidth();
	iHeight = pInputImage->GetHeight();

    listPaths.clear();

    GeoComp::Grid8<int> grid(iHeight, iWidth);
    GeoComp::FM2<float, GeoComp::Grid8<int> > fastMarching(grid);

    // Combined action map and its first-order derivatives
    CArray2D<float> arrayAction, arrayActionDerX, arrayActionDerY;

    CArray2D<int> arrayVoronoi;
    float *U, *DUx, *DUy, *P;
    int *V;
    int sources_index[2];

    sources_index[0] = grid.index((int)(pfStart.y + 0.5f), (int)(pfStart.x + 0.5f));
    sources_index[1] = grid.index((int)(pfEnd.y + 0.5f), (int)(pfEnd.x + 0.5f));

    arrayAction.Init(iWidth, iHeight);
    arrayActionDerX.Init(iWidth, iHeight);
    arrayActionDerY.Init(iWidth, iHeight);
    arrayVoronoi.Init(iWidth, iHeight);

    U = arrayAction.GetBuffer();
    DUx = arrayActionDerX.GetBuffer();
    DUy = arrayActionDerY.GetBuffer();
    V = arrayVoronoi.GetBuffer();

    P = arrayPotential.GetBuffer();

    fastMarching.voronoiFast(sources_index, 2, P, V, U, DUy, DUx);

    arrayAction.FillBoundary(numeric_limits<float>::max());

    #ifdef GRAPHICWND_H
    static CArray2DFloatWnd *pWndU = NULL;
    if (pWndU==NULL)
        pWndU = new CArray2DFloatWnd(iWidth, iHeight, "U");
    pWndU->SetArray(arrayAction);
    pWndU->UpdateImage();
    pWndU->show();

    static CArray2DFloatWnd *pWndVoronoi = NULL;
    CArray2D<float> arrayVoronoiFlt;
    arrayVoronoiFlt.InitCast(arrayVoronoi);
    if (pWndVoronoi==NULL)
        pWndVoronoi = new CArray2DFloatWnd(iWidth, iHeight, "Voronoi");
    pWndVoronoi->SetArray(arrayVoronoiFlt);
    pWndVoronoi->UpdateImage();
    pWndVoronoi->show();
    #endif

    //list<MySaddlePair> listSaddlePairs;
    //MySaddlePair saddlePair;

    voronoiGraph.SetActionAndVoronoiMaps(&arrayAction, &arrayVoronoi);

    // Get the point located on the Voronoi edge having the lowest action
    CCouple<int> vorPtStart;

    if (voronoiGraph.GetVoronoiEdgePointMinimalAction(vorPtStart)==false)
    {
        cout<<"WARNING in CPiecewiseGeodesicCombinationBase::ComputeAdmissiblePaths(...): ";
        cout<<" No suitable Voronoi edge point found. List of admissible paths will be empty"<<endl;
        return;
    }

    #ifdef GRAPHICWND_H
    CArray2D<float> arrayVoronoiEdges;
    arrayVoronoiEdges.InitCast(voronoiGraph.arrayPointelsOnVoronoiEdges);

    static CArray2DFloatWnd *pWndVoronoiEdges = NULL;
    if (pWndVoronoiEdges==NULL)
        pWndVoronoiEdges = new CArray2DFloatWnd(iWidth, iHeight, "Voronoi edges");
    pWndVoronoiEdges->SetArray(arrayVoronoiEdges);
    pWndVoronoiEdges->UpdateImage();
    pWndVoronoiEdges->show();
    #endif

    // Compute digital 4-connected Voronoi medial curve and smooth action values along this curve
    voronoiGraph.GetEdgeFromPointel(vorPtStart, pathVoronoi);
    for (int iSmooth=0; iSmooth<10; iSmooth++)
        pathVoronoi.SmoothActionsAlongEdge();

    // Get local minima
    list<CVoronoiEdgePoint> listSaddles;
    list<CVoronoiEdgePoint>::const_iterator itSaddle;

    pathVoronoi.ComputeLocalMinimaAlongEdge(listSaddles);
    if (listSaddles.size()==0)
    {
        cout<<"WARNING in CPiecewiseGeodesicCombinationBase::ComputeAdmissiblePaths(...): ";
        cout<<" No saddle point found. List of admissible paths will be empty"<<endl;
        return;
    }

    CCouple<float> pfCurrent, pfLeft, pfRight, pfSaddleStart, pfSaddleEnd, pfHalf(0.5f, 0.5f);
    int iCellLeft, iCellRight;
    unsigned int iNbPathsBuilt = 0;
    CRealPath2D pathTemp, pathNew;
    bool bBuildPath;

    for (itSaddle=listSaddles.begin(); itSaddle!=listSaddles.end() && iNbPathsBuilt<iNbAdmissiblePathsMax; itSaddle++)
    {
        // The saddle point is located on a pointel
        // So, convert pointel coordinates into real pixel coordinates
        pfCurrent = (CCouple<float>)(itSaddle->piPos) - pfHalf;

        pfLeft = pfCurrent + itSaddle->vfNormal;
        pfRight = pfCurrent - itSaddle->vfNormal;

        iCellLeft = arrayVoronoi.Element((CCouple<int>)(pfLeft + pfHalf));
        iCellRight = arrayVoronoi.Element((CCouple<int>)(pfRight + pfHalf));

        if (iCellLeft!=iCellRight)
        {
            if (iCellLeft<iCellRight)
            {
                pfSaddleStart = pfLeft;
                pfSaddleEnd = pfRight;
            }
            else {
                pfSaddleStart = pfRight;
                pfSaddleEnd = pfLeft;
            }

            bBuildPath = GradientDescentActionMap(pfSaddleStart, arrayActionDerX, arrayActionDerY, pathTemp);
            if (bBuildPath==true)
            {
                pathTemp.reverse();
                if ((pathTemp.front()-pfStart).L2Norm2()>1.0f)
                {
                    cout<<"Start point was not reached by gradient descent"<<endl;
                    bBuildPath = false;
                }
            }

            if (bBuildPath==true)
            {
                pathNew = pathTemp;
                bBuildPath = GradientDescentActionMap(pfSaddleEnd, arrayActionDerX, arrayActionDerY, pathTemp);
            }

            if (bBuildPath==true)
            {
                if ((pathTemp.back()-pfEnd).L2Norm2()>1.0f)
                {
                    cout<<"End point was not reached by gradient descent"<<endl;
                    bBuildPath = false;
                }
            }

            if (bBuildPath==true)
            {
                // Assemble the two geodesic curves and add the resulting
                // piecewise-geodesic curve to the list of admissible paths
                pathNew.splice(pathNew.end(), pathTemp);
                listPaths.push_back(pathNew);
                iNbPathsBuilt++;
            }
        }
        else {
            cout<<"Wrong saddle point. Points on the left and right of the saddle point are in the same Voronoi cell"<<endl;
        }
    }
}

void CPiecewiseGeodesicCombinationBase::FindBestCombinationGreedy()
{
	list<list<CRealPath2D> >::const_iterator itListPaths;
	list<CRealPath2D>::const_iterator itPath;
	vector<unsigned int> vectCurrentCombination, vectBestCombinationIter, vectCandidateCombination, listNbPaths;
	unsigned int iVertex, iNbVertices;
	bool bStop, bIterWrongEnergies;
	float fEnergy, fEnergyMin, fEnergyMinIter;

	iNbVertices = listVertices.size();

	vectCurrentCombination.resize(iNbVertices);
    vectBestCombination.resize(iNbVertices);
	vectBestCombinationIter.resize(iNbVertices);
	vectCandidateCombination.resize(iNbVertices);

    listNbPaths.resize(iNbVertices); // Number of paths per admissible set

	for (iVertex=0; iVertex<iNbVertices; iVertex++)
	{
        vectCurrentCombination[iVertex] = 0;
        vectBestCombination[iVertex] = 0;
        vectBestCombinationIter[iVertex] = 0;
    }
	bStop = false;

    for (iVertex=0, itListPaths=listListsAdmissiblePaths.begin(); iVertex<iNbVertices; iVertex++, itListPaths++)
        listNbPaths[iVertex] = itListPaths->size();

    fEnergyMin = GetEnergy(vectCurrentCombination);

    cout<<"GREEDY LOCAL SEARCH"<<endl<<endl;

	while (bStop==false)
	{
	    bStop = true;
	    vectCurrentCombination = vectBestCombinationIter;

        fEnergyMinIter = numeric_limits<float>::max();
        bIterWrongEnergies = true;

	    // Loop over candidate combinations
	    for (iVertex=0; iVertex<iNbVertices; iVertex++)
	    {
            vectCandidateCombination = vectCurrentCombination;
            if (vectCandidateCombination[iVertex]<listNbPaths[iVertex]-1)
            {
                bStop = false;

                vectCandidateCombination[iVertex]++;
                fEnergy = GetEnergy(vectCandidateCombination);

                if (fEnergy<fEnergyMinIter)
                {
                    bIterWrongEnergies = false;

                    fEnergyMinIter = fEnergy;
                    vectBestCombinationIter = vectCandidateCombination;

                    if (fEnergy<fEnergyMin)
                    {
                        fEnergyMin = fEnergy;
                        vectBestCombination = vectCandidateCombination;

                        cout<<"Found better combination ";

                        for (iVertex=0; iVertex<iNbVertices; iVertex++)
                            cout<<vectBestCombination[iVertex]<<",";
                        cout<<" E="<<fEnergyMin;

                        if (fWeightSimplicity!=0.0f)
                            cout<<"  Esimplicity="<<fEnergySimplicity;
                        if (fWeightEdge!=0.0f)
                            cout<<"  Eedge="<<fEnergyEdge;
                        if (fWeightRegion!=0.0f)
                            cout<<"  Eregion="<<fEnergyRegion;
                        cout<<endl;
                    }
                }
            }
	    }

	    if (bIterWrongEnergies==true)
	    {
            float fInnerRegionArea, fInnerRegionAreaMin;

            cout<<"All candidate combinations have infinite energy. Select combination with minimal area from "; //<<endl;
            for (iVertex=0; iVertex<iNbVertices; iVertex++)
                cout<<vectCurrentCombination[iVertex]<<",";
            cout<<endl;

            fInnerRegionAreaMin = numeric_limits<float>::max();

            // Loop over candidate combinations
            for (iVertex=0; iVertex<iNbVertices; iVertex++)
            {
                vectCandidateCombination = vectCurrentCombination;
                if (vectCandidateCombination[iVertex]<listNbPaths[iVertex]-1)
                {
                    vectCandidateCombination[iVertex]++;
                    fInnerRegionArea = GetSignedArea(vectCandidateCombination);

                    if (fInnerRegionArea<fInnerRegionAreaMin)
                    {
                        fInnerRegionAreaMin = fInnerRegionArea;
                        vectBestCombinationIter = vectCandidateCombination;
                    }
                }
            }

            cout<<"Combination of minimal area: ";
            for (iVertex=0; iVertex<iNbVertices; iVertex++)
                cout<<vectBestCombinationIter[iVertex]<<",";
            cout<<endl;
	    }
	}
}

void CPiecewiseGeodesicCombinationBase::FindBestCombinationExhaustive()
{
	list<list<CRealPath2D> >::const_iterator itListPaths;
	vector<unsigned int> vectCurrentCombination;
	unsigned int iVertex, iNbVertices;
	bool bAllIndicesTested;
	unsigned int iCurrentIndex;
	float fEnergy, fEnergyMin;

	iNbVertices = listVertices.size();

	vectCurrentCombination.resize(iNbVertices);
	vectBestCombination.resize(iNbVertices);

	fEnergyMin = numeric_limits<float>::max();

	for (iVertex=0; iVertex<iNbVertices; iVertex++)
        vectCurrentCombination[iVertex] = 0;

    cout<<"EXHAUSTIVE SEARCH"<<endl<<endl;

	bAllIndicesTested = false;

	while (bAllIndicesTested==false)
	{
		fEnergy = GetEnergy(vectCurrentCombination);

		if (fEnergy<fEnergyMin)
		{
			fEnergyMin = fEnergy;
			vectBestCombination = vectCurrentCombination;

			cout<<"Found better combination ";

			for (iVertex=0; iVertex<iNbVertices; iVertex++)
				cout<<vectBestCombination[iVertex]<<",";
			cout<<" E="<<fEnergyMin;

            if (fWeightSimplicity!=0.0f)
                cout<<"  Esimplicity="<<fEnergySimplicity;
            if (fWeightEdge!=0.0f)
                cout<<"  Eedge="<<fEnergyEdge;
            if (fWeightRegion!=0.0f)
                cout<<"  Eregion="<<fEnergyRegion;
            cout<<endl;
		}

		iCurrentIndex = 0;
		vectCurrentCombination[iCurrentIndex]++;
        itListPaths = listListsAdmissiblePaths.begin();

		while (iCurrentIndex<iNbVertices && vectCurrentCombination[iCurrentIndex]==itListPaths->size())
		{
			vectCurrentCombination[iCurrentIndex] = 0;
			iCurrentIndex++;
			itListPaths++;
			if (iCurrentIndex<iNbVertices)
				vectCurrentCombination[iCurrentIndex]++;
		}
		if (iCurrentIndex==iNbVertices)
			bAllIndicesTested = true;
	}
}

void CPiecewiseGeodesicCombinationBase::FindBestCombinationTwoVertices()
{
	/*
	list<const CRealPath2D *> listPathsSetToTest;
	list<CRealPath2D> *pListAd;
	list<CRealPath2D>::const_iterator itPath, itPathOther, itPathBest, itPathOtherBest;
	unsigned int iPath, iPathOther;

	// int iIter, iSmooth;
	float fEnergy, fEnergyMin;

	fEnergyMin = numeric_limits<float>::max();
	listPathsSetToTest.push_back(NULL);
	listPathsSetToTest.push_back(NULL);
    pListAd = &(listListsAdmissiblePaths.front());

	for (itPath=pListAd->begin(), iPath=0; iPath<pListAd->size()-1; itPath++, iPath++)
	{
		// listIndices[0] = iPath;
		listPathsSetToTest.front() = &(*itPath);

        itPathOther = itPath;
		for (itPathOther++, iPathOther=iPath+1; iPathOther<pListAd->size(); itPathOther++, iPathOther++)
		{
			// listIndices[1] = iPathOther;
			listPathsSetToTest.back() = &(*itPathOther);

			cout<<"Config "<<iPath<<","<<iPathOther<<endl;
			fEnergy = GetEnergy(&listPathsSetToTest);

			if (fEnergy<fEnergyMin)
			{
				fEnergyMin = fEnergy;
				itPathBest = itPath;
				itPathOtherBest = itPathOther;
				cout<<"Best config E="<<fEnergyMin<<endl;
			}
		}
	}

	listSelectedPaths.clear();
	listSelectedPaths.front() = *itPathBest;
	listSelectedPaths.back() = *itPathOtherBest;
	*/

	const list<CRealPath2D> *pListPaths1, *pListPaths2;
	vector<unsigned int> vectCurrentCombination;
	unsigned int iNbPaths;
	float fEnergy, fEnergyMin;

	vectCurrentCombination.resize(2);
	vectBestCombination.clear();

    pListPaths1 = &(listListsAdmissiblePaths.front());
    pListPaths2 = &(listListsAdmissiblePaths.back());

    if (pListPaths1->size()!=pListPaths2->size())
    {
        cerr<<"ERROR : the two lists of admissible paths have different sizes"<<endl;
        return;
    }
    iNbPaths = pListPaths1->size();

	fEnergyMin = numeric_limits<float>::max();

    cout<<"SEARCH WITH TWO VERTICES"<<endl<<endl;

	for (vectCurrentCombination[0]=1; vectCurrentCombination[0]<iNbPaths; vectCurrentCombination[0]++)
	{
        for (vectCurrentCombination[1]=iNbPaths-vectCurrentCombination[0]; vectCurrentCombination[1]<iNbPaths; vectCurrentCombination[1]++)
        {
            fEnergy = GetEnergy(vectCurrentCombination);
            if (fEnergy<fEnergyMin)
            {
                fEnergyMin = fEnergy;
                vectBestCombination = vectCurrentCombination;

                cout<<"Found better combination "<<vectBestCombination[0]<<","<<vectBestCombination[1];
                cout<<" E="<<fEnergyMin;

                if (fWeightSimplicity!=0.0f)
                    cout<<"  Esimplicity="<<fEnergySimplicity;
                if (fWeightEdge!=0.0f)
                    cout<<"  Eedge="<<fEnergyEdge;
                if (fWeightRegion!=0.0f)
                    cout<<"  Eregion="<<fEnergyRegion;
                cout<<endl;
            }
        }
    }
}

// Puts the innermost first, outermost last ?
void CPiecewiseGeodesicCombinationBase::SortAdmissiblePathsByExteriority()
{
    list<list<CRealPath2D> >::iterator itListAdPaths;
    list<CRealPath2D>::iterator itPath;
    list<CRealPath2D> *pList;

    list<list<CIntegerPath2D> >::iterator itListAdPathsDigital;
    list<CIntegerPath2D>::iterator itPathDigital;
    list<CIntegerPath2D> *pListDigital;

    CCouple<float> pfStart, pfEnd;

    if (listListsAdmissiblePaths.size()!=listListsAdmissiblePathsDigital.size())
    {
        cerr<<"ERROR in CPiecewiseGeodesicCombinationBase::SortAdmissiblePathsByExteriority():";
        cerr<<" lists of real and digital admissible sets have different sizes"<<endl;
        return;
    }

    for (itListAdPaths=listListsAdmissiblePaths.begin(), itListAdPathsDigital=listListsAdmissiblePathsDigital.begin();
        itListAdPaths!=listListsAdmissiblePaths.end(); itListAdPaths++, itListAdPathsDigital++)
    {
        pList = &(*itListAdPaths);
        pListDigital = &(*itListAdPathsDigital);

        if (pList->size()!=pListDigital->size())
        {
            cerr<<"ERROR in CPiecewiseGeodesicCombinationBase::SortAdmissiblePathsByExteriority():";
            cerr<<" real and digital admissible sets have different sizes"<<endl;
            return;
        }

        for (itPath=pList->begin(), itPathDigital=pListDigital->begin(); itPath!=pList->end(); itPath++, itPathDigital++)
        {
            itPath->fSortingValue = itPath->GetSignedAreaWithLineSegment();
            itPathDigital->fSortingValue = itPath->fSortingValue;
        }
        pList->sort();
        pListDigital->sort();

        CImage2DFloatRGBPixel fRgbStart(0,0,255), fRgbEnd(255,0,0);
        float fCoef;
        int iPath;

        for (iPath=0, itPath=pList->begin(), itPathDigital=pListDigital->begin(); itPath!=pList->end(); iPath++, itPath++, itPathDigital++)
        {
            fCoef = (float)iPath/(float)(pList->size()-1);
            itPath->rgbEdges = (fCoef*fRgbEnd + (1.0f-fCoef)*fRgbStart).ToByteRGBPixel();
            itPathDigital->rgbEdges = (fCoef*fRgbEnd + (1.0f-fCoef)*fRgbStart).ToByteRGBPixel();
        }
    }
}

/*
void CPiecewiseGeodesicCombinationBase::SelectPathsByIndex(const CArray1D<unsigned int> &arrayIndices)
{
	list<list<CRealPath2D> >::const_iterator itListPaths;
	list<CRealPath2D>::const_iterator itPath;
	list<list<CIntegerPath2D> >::const_iterator itListDigitalPaths;
	list<CIntegerPath2D>::const_iterator itDigitalPath;
	unsigned int iVertex, iNbVertices;
	float fEnergy;

    iNbVertices = listListsAdmissiblePaths.size();

    listChosenPaths.clear();
    listChosenPathsDigital.clear();
    for (iVertex=0, itListPaths=listListsAdmissiblePaths.begin(), itListDigitalPaths=listListsAdmissiblePathsDigital.begin();
    iVertex<iNbVertices; iVertex++, itListPaths++, itListDigitalPaths++)
    {
        unsigned int iInd;
        for (iInd=0, itPath=itListPaths->begin(), itDigitalPath=itListDigitalPaths->begin();
            iInd<arrayIndices[iVertex]; iInd++)
        {
            itPath++;
            itDigitalPath++;
        }
        listChosenPaths.push_back(*itPath);
        listChosenPathsDigital.push_back(*itDigitalPath);
    }
    arrayChosenPathsIndices = arrayIndices;

    fEnergy = GetEnergyByIndex(CArray1D<unsigned int>());
    cout<<"E="<<fEnergy;

    if (fWeightEdge!=0.0f)
        cout<<"  Egradient="<<fEnergyEdge;
    if (fWeightSimplicity!=0.0f)
        cout<<"  Eoverlap="<<fEnergySimplicity;
    if (fWeightRegion!=0.0f)
        cout<<"  Eregion="<<fEnergyRegion;
    cout<<endl;
}*/

float CPiecewiseGeodesicCombinationBase::Overlap(const CIntegerPath2D &intPath1, const CIntegerPath2D &intPath2) const
{
    CIntegerPath2D::const_iterator itPoint;
	CArray2D<float> arrayWeight1; //, arrayWeight2;
	float fOverlap2To1; //, fOverlap1To2;
	float *pWeight;
	int iWidth; //, iPoint;

	// Overlapping "from path 2 to path 1"
	arrayWeight1.Init(pInputImage->GetSize());
	arrayWeight1.Fill(0.0f);
	iWidth = arrayWeight1.GetWidth();

	for (itPoint=intPath1.begin(); itPoint!=intPath1.end(); itPoint++)
	{
        pWeight = arrayWeight1.GetBuffer() + arrayWeight1.GetOffset(*itPoint);
        pWeight[0] = 1.0f;
        pWeight[-1] = 1.0f;
        pWeight[1] = 1.0f;
        pWeight[-iWidth] = 1.0f;
        pWeight[iWidth] = 1.0f;
	}

	fOverlap2To1 = 0.0f;
	for (itPoint=intPath2.begin(); itPoint!=intPath2.end(); itPoint++)
		fOverlap2To1 += arrayWeight1.Element(*itPoint);

	return fOverlap2To1;
}

class CPointAndIntersection
{
  public:
    CCouple<int> piPos;
    int iTarget; // Index of point crossed, -1 if the point is of multiplicity 1
    enum {INTERSECTING, INTERSECTED} role;
    enum {ENTER, LEAVE, ENTER_AND_LEAVE} type;
    enum {RIGHT, LEFT, HALFTURN} side;
  public:
    CPointAndIntersection() {}
};

float CPiecewiseGeodesicCombinationBase::InvertedAreaClosed(const CIntegerPath2D &intPath) const
{
    CIntegerPath2D::const_iterator itPoint;
    CArray2D<unsigned int> arrayMultiplicity;
    CArray2D<int> arrayIndices;
    CCouple<int> viNormal, viCrossedNormal, viTangent;;
    unsigned int iPoint, iStart, iNbPoints;

    // Curve with intersection data
    CPointAndIntersection ptAndInter;
    vector<CPointAndIntersection> path;

    // Total inverted area, as defined in Eq. (11)
    float fInvertedArea;

    // Area of the last found single or double loop
    float fArea;

    arrayMultiplicity.Init(pInputImage->GetSize());
    arrayMultiplicity.Fill(0);

    arrayIndices.Init(pInputImage->GetSize());
    arrayIndices.Fill(-1);

    iNbPoints = intPath.size();
    path.reserve(iNbPoints);

    unsigned int iPointInit[4]; // Left, right, top and bottom
    unsigned int iInit = 0;

    for (iInit=0; iInit<4; iInit++)
        iPointInit[iInit] = 0;

    for (itPoint=intPath.begin(), iPoint=0; iPoint<iNbPoints; itPoint++, iPoint++)
    {
        ptAndInter.piPos = *itPoint;
        ptAndInter.iTarget = -1;

        path.push_back(ptAndInter);

        arrayMultiplicity.Element(*itPoint)++;

        // If a point with multiplicity>2 is met, twisting measure is infinite
        if (arrayMultiplicity.Element(*itPoint)>2)
            return numeric_limits<float>::max();

        if (path.back().piPos.x < path[iPointInit[0]].piPos.x)
            iPointInit[0] = iPoint;
        else if (path.back().piPos.x > path[iPointInit[1]].piPos.x)
            iPointInit[1] = iPoint;

        if (path.back().piPos.y < path[iPointInit[2]].piPos.y)
            iPointInit[2] = iPoint;
        else if (path.back().piPos.y > path[iPointInit[3]].piPos.y)
            iPointInit[3] = iPoint;
    }

    bool bStartPointFound = false;

    // Find a point of multiplicity 1 with normal vector pointing inward
    iInit = 0;
    while (iInit<4 && bStartPointFound==false)
    {
        iPoint = iPointInit[iInit];

        // pPtAndInter = &(path[iPointInit[iInit]]);
        if (arrayMultiplicity.Element(path[iPoint].piPos)==1)
        {
            viNormal = (path[(iPoint+1)%iNbPoints].piPos - path[(iPoint+iNbPoints-1)%iNbPoints].piPos).Perp();
            if ((iInit==0 && viNormal.x>0) || (iInit==1 && viNormal.x<0)
                || (iInit==2 && viNormal.y>0) || (iInit==3 && viNormal.y<0))
            {
                iStart = iPoint;
                bStartPointFound = true;
            }
        }
        iInit++;
    }

    if (bStartPointFound==false)
    {
        CCouple<int> piCurrent;

        // Try other initial points
        arrayIndices.Fill(-1);
        for (iPoint = 0; iPoint<path.size(); iPoint++)
        {
            if (arrayIndices.Element(path[iPoint].piPos)==-1)
                arrayIndices.Element(path[iPoint].piPos) = (int)iPoint;
        }

        // From top
        piCurrent.x = (path[iPointInit[0]].piPos.x + path[iPointInit[1]].piPos.x)/2;
        piCurrent.y = 0;
        while (piCurrent.y<pInputImage->GetHeight() && arrayMultiplicity.Element(piCurrent)==0)
            piCurrent.y++;
        if (piCurrent.y<pInputImage->GetHeight() && arrayMultiplicity.Element(piCurrent)==1)
        {
            iPoint = (unsigned int)arrayIndices.Element(piCurrent);
            viNormal = (path[(iPoint+1)%iNbPoints].piPos - path[(iPoint+iNbPoints-1)%iNbPoints].piPos).Perp();
            if (viNormal.y>0)
            {
                iStart = iPoint;
                bStartPointFound = true;
            }
        }

        // From bottom
        if (bStartPointFound==false)
        {
            piCurrent.x = (path[iPointInit[0]].piPos.x + path[iPointInit[1]].piPos.x)/2;
            piCurrent.y = pInputImage->GetHeight()-1;
            while (piCurrent.y>=0 && arrayMultiplicity.Element(piCurrent)==0)
                piCurrent.y--;

            if (piCurrent.y>=0 && arrayMultiplicity.Element(piCurrent)==1)
            {
                iPoint = (unsigned int)arrayIndices.Element(piCurrent);
                viNormal = (path[(iPoint+1)%iNbPoints].piPos - path[(iPoint+iNbPoints-1)%iNbPoints].piPos).Perp();
                if (viNormal.y<0)
                {
                    iStart = iPoint;
                    bStartPointFound = true;
                }
            }
        }
    }

    if (bStartPointFound==false)
    {
        cout<<"WARNING in CPiecewiseGeodesicCombinationBase::InvertedAreaClosed(...): ";
        cout<<" failed to find a suitable starting point. Returning infinity..."<<endl;
        return numeric_limits<float>::max();
    }

    int iDotProduct;
    unsigned int iPointIntersected, iPass, iPointPred;;
    bool bOnCurve;

    fInvertedArea = 0.0f;
    iPoint = iStart;
    bOnCurve = false;

    for (iPass=0; iPass<iNbPoints; iPoint++, iPass++)
    {
        if (iPoint==iNbPoints)
            iPoint = 0;

        if (bOnCurve==false)
        {
            if (arrayIndices.Element(path[iPoint].piPos)==-1)
                arrayIndices.Element(path[iPoint].piPos) = (int)iPoint;
            else {
                // Entering on curve, check side...
                bOnCurve = true;

                iPointIntersected = (unsigned int)arrayIndices.Element(path[iPoint].piPos);
                path[iPoint].iTarget = (int)iPointIntersected;
                path[iPoint].role = CPointAndIntersection::INTERSECTING;
                path[iPoint].type = CPointAndIntersection::ENTER;

                path[iPointIntersected].iTarget = (int)iPoint;
                path[iPointIntersected].role = CPointAndIntersection::INTERSECTED;
                path[iPointIntersected].type = CPointAndIntersection::ENTER;

                viCrossedNormal = (path[(iPointIntersected+1)%iNbPoints].piPos -
                            path[(iPointIntersected+iNbPoints-1)%iNbPoints].piPos).Perp();

                viTangent = path[iPoint].piPos - path[(iPoint+iNbPoints-1)%iNbPoints].piPos;

                // Let u and v be the arc lengths of the crossed point and current point, respectively
                // If n(u).t(v)>0, C(v) enters on C(u) from the left
                // See Section 4.2
                iDotProduct = DotProduct(viCrossedNormal, viTangent);
                if (iDotProduct>0)
                {
                    path[iPoint].side = CPointAndIntersection::LEFT;
                    path[iPointIntersected].side = CPointAndIntersection::LEFT;
                }
                else if (iDotProduct<0)
                {
                    path[iPoint].side = CPointAndIntersection::RIGHT;
                    path[iPointIntersected].side = CPointAndIntersection::RIGHT;
                }
                else {
                    path[iPoint].side = CPointAndIntersection::HALFTURN;
                    path[iPointIntersected].side = CPointAndIntersection::HALFTURN;
                }

                if (path[iPoint].side==CPointAndIntersection::LEFT)
                {
                    // Go along curve from the intersected point to the intersecting point
                    CIntegerPath2D invPath;
                    unsigned int invPoint;
                    bool bLoopFound = false;

                    invPath.push_back(path[iPointIntersected].piPos);
                    invPoint = iPointIntersected+1;
                    if (invPoint==iNbPoints)
                        invPoint = 0;
                    while (bLoopFound==false)
                    {
                        if (invPoint==iPoint || invPoint==iPointIntersected)
                            bLoopFound = true; // Perhaps a single loop?
                        else if (path[invPoint].iTarget!=-1
                            && (path[invPoint].role==CPointAndIntersection::INTERSECTING || path[invPoint].role==CPointAndIntersection::INTERSECTED)
                            && path[invPoint].side==CPointAndIntersection::LEFT)
                            bLoopFound = true; // Perhaps a double loop?
                        else
                        {
                            invPath.push_back(path[invPoint].piPos);
                            invPoint++;
                            if (invPoint==iNbPoints)
                                invPoint = 0;
                        }
                    }

                    if (invPoint==iPoint)
                    {
                        // Single loop found
                        fArea = invPath.GetArea4connected();
                        if (fArea<-10.0f)
                        {
                            fInvertedArea += fArea;
                            // invPath.rgbEdges.Set(0,0,255);
                            // invPath.DrawInImage(imgDraw, iZoom);
                        }
                    }
                    else if (invPoint==iPointIntersected)
                    {
                        cout<<"WARNING in CPiecewiseGeodesicCombinationBase::InvertedAreaClosed(...): ";
                        cout<<" gone over the entire curve without detecting loop"<<endl;
                    }
                    else if (path[invPoint].iTarget!=-1
                            && (path[invPoint].role==CPointAndIntersection::INTERSECTING || path[invPoint].role==CPointAndIntersection::INTERSECTED)
                            && path[invPoint].side==CPointAndIntersection::LEFT)
                    {
                        if (path[invPoint].role==CPointAndIntersection::INTERSECTED)
                        {
                            invPoint = (unsigned int)(path[invPoint].iTarget);

                            invPath.push_back(path[invPoint].piPos);
                            if (++invPoint==iNbPoints)
                                invPoint = 0;
                            bLoopFound = false;

                            while (bLoopFound==false)
                            {
                                if (invPoint==iPoint || invPoint==iPointIntersected)
                                    bLoopFound = true; // Perhaps a double loop?
                                else if (path[invPoint].iTarget!=-1
                                    && (path[invPoint].role==CPointAndIntersection::INTERSECTING || path[invPoint].role==CPointAndIntersection::INTERSECTED)
                                    && path[invPoint].side==CPointAndIntersection::LEFT)
                                    bLoopFound = true; // This should normally not happen
                                else
                                {
                                    invPath.push_back(path[invPoint].piPos);
                                    invPoint++;
                                    if (invPoint==iNbPoints)
                                        invPoint = 0;
                                }
                            }
                            if (invPoint==iPoint)
                            {
                                // Double loop found
                                fArea = invPath.GetArea4connected();
                                if (fArea<-10.0f)
                                {
                                    fInvertedArea += fArea;
                                    // invPath.rgbEdges.Set(0,0,255);
                                    // invPath.DrawInImage(imgDraw, iZoom);
                                }
                            }
                            else {
                                cout<<"WARNING in CPiecewiseGeodesicCombinationBase::InvertedAreaClosed(...): ";
                                cout<<" gone over the entire curve without detecting loop"<<endl;
                            }
                        }
                    }
                }
            }
        }
        else {
            if (arrayIndices.Element(path[iPoint].piPos)==-1)
            {
                // Leaving curve, check side...
                arrayIndices.Element(path[iPoint].piPos) = (int)iPoint;
                bOnCurve = false;

                iPointPred = (iPoint+iNbPoints-1)%iNbPoints;
                iPointIntersected = (unsigned int)arrayIndices.Element(path[iPointPred].piPos);

                if (path[iPointPred].iTarget==-1)
                {
                    path[iPointPred].iTarget = (int)iPointIntersected;
                    path[iPointPred].role = CPointAndIntersection::INTERSECTING;
                    path[iPointPred].type = CPointAndIntersection::LEAVE;

                    path[iPointIntersected].iTarget = (int)iPointPred;
                    path[iPointIntersected].role = CPointAndIntersection::INTERSECTED;
                    path[iPointIntersected].type = CPointAndIntersection::LEAVE;

                    viCrossedNormal = (path[(iPointIntersected+1)%iNbPoints].piPos -
                                path[(iPointIntersected+iNbPoints-1)%iNbPoints].piPos).Perp();

                    viTangent = path[iPoint].piPos - path[iPointPred].piPos;

                    // Let u and v be the arc lengths of the crossed point and current point, respectively
                    // If n(u).t(v)<0, C(v) leaves from C(u) to go left
                    // See Section 4.2
                    iDotProduct = DotProduct(viCrossedNormal, viTangent);
                    if (iDotProduct<0)
                    {
                        path[iPointPred].side = CPointAndIntersection::LEFT;
                        path[iPointIntersected].side = CPointAndIntersection::LEFT;
                    }
                    else if (iDotProduct>0)
                    {
                        path[iPointPred].side = CPointAndIntersection::RIGHT;
                        path[iPointIntersected].side = CPointAndIntersection::RIGHT;
                    }
                    else {
                        path[iPointPred].side = CPointAndIntersection::HALFTURN;
                        path[iPointIntersected].side = CPointAndIntersection::HALFTURN;
                    }
                }
                else {
                    path[iPointPred].type = CPointAndIntersection::ENTER_AND_LEAVE;
                    path[iPointIntersected].type = CPointAndIntersection::ENTER_AND_LEAVE;
                }
            }
        }
    }

    return fabs(fInvertedArea);
}

void CPiecewiseGeodesicCombinationBase::ConvertRealPathToDigitalPathSimple(const CRealPath2D &realPath, CIntegerPath2D &intPath) const
{
    CRealPath2D::const_iterator itRealPoint;
    const CCouple<int> *pPoint;
    CArray2D<CCouple<int> *> arrayPtrPoints;
    CCouple<float> pfHalf(0.5f, 0.5f);
    CCouple<int> piLastAdded, piCurrent;
    CCouple<int> viAbsDiff;

    arrayPtrPoints.Init(pInputImage->GetWidth(), pInputImage->GetHeight());
    arrayPtrPoints.Fill(NULL);

    intPath.clear();

    itRealPoint = realPath.begin();
    piLastAdded = (CCouple<int>)(*itRealPoint + pfHalf);
    if (arrayPtrPoints.Element(piLastAdded)==NULL)
    {
        intPath.push_back(piLastAdded);
        arrayPtrPoints.Element(piLastAdded) = &intPath.back();
    }

    for (itRealPoint++; itRealPoint!=realPath.end(); itRealPoint++)
    {
        piCurrent = (CCouple<int>)(*itRealPoint + pfHalf);
        if (piCurrent!=piLastAdded)
        {
            viAbsDiff.x = abs(piCurrent.x-piLastAdded.x);
            viAbsDiff.y = abs(piCurrent.y-piLastAdded.y);

            if (viAbsDiff.x+viAbsDiff.y==1)
            {
                if (arrayPtrPoints.Element(piCurrent)==NULL)
                {
                    intPath.push_back(piCurrent);
                    arrayPtrPoints.Element(piCurrent) = &intPath.back();
                    piLastAdded = piCurrent;
                }
                else {
                    pPoint = arrayPtrPoints.Element(piCurrent);
                    while (pPoint!=&intPath.back())
                    {
                        arrayPtrPoints.Element(intPath.back()) = NULL;
                        intPath.pop_back();
                    }
                    piLastAdded = piCurrent;
                }
            }
            else {
                // Draw digital line from piLastAdded to piLastCurrent
                CIntegerPath2D pathDigitalLineSeg;
                CIntegerPath2D::const_iterator itPointSeg;

                if (viAbsDiff.x==1 && viAbsDiff.y==1)
                {
                    pathDigitalLineSeg.push_back(piLastAdded);
                    pathDigitalLineSeg.push_back(CCouple<int>(piLastAdded.x, piCurrent.y));
                    pathDigitalLineSeg.push_back(piCurrent);
                }
                else {
                    pathDigitalLineSeg.SetDigitalLineSegment4connected(piLastAdded, piCurrent);
                }

                for (itPointSeg=++pathDigitalLineSeg.begin(); itPointSeg!=pathDigitalLineSeg.end(); itPointSeg++)
                {
                    piCurrent = *itPointSeg;
                    if (arrayPtrPoints.Element(piCurrent)==NULL)
                    {
                        intPath.push_back(piCurrent);
                        arrayPtrPoints.Element(piCurrent) = &intPath.back();
                        piLastAdded = piCurrent;
                    }
                    else {
                        pPoint = arrayPtrPoints.Element(piCurrent);
                        while (pPoint!=&intPath.back())
                        {
                            arrayPtrPoints.Element(intPath.back()) = NULL;
                            intPath.pop_back();
                        }
                        piLastAdded = piCurrent;
                    }
                }
            }
        }
    }

    intPath.rgbEdges = realPath.rgbEdges;
    intPath.rgbVertices = realPath.rgbVertices;
    intPath.bDisplayEdges = realPath.bDisplayVertices;
    intPath.bDisplayEdges = realPath.bDisplayEdges;
}

void CPiecewiseGeodesicCombinationBase::ConvertRealPathToDigitalPath(const CRealPath2D &realPath, CIntegerPath2D &intPath) const
{
    CRealPath2D::const_iterator itRealPoint;
    CCouple<float> pfHalf(0.5f, 0.5f);
    CCouple<int> piLastAdded, piCurrent;
    CCouple<int> viAbsDiff;

    intPath.clear();

    itRealPoint = realPath.begin();
    piLastAdded = (CCouple<int>)(*itRealPoint + pfHalf);
    intPath.push_back(piLastAdded);

    for (itRealPoint++; itRealPoint!=realPath.end(); itRealPoint++)
    {
        piCurrent = (CCouple<int>)(*itRealPoint + pfHalf);
        if (piCurrent!=piLastAdded)
        {
            viAbsDiff.x = abs(piCurrent.x-piLastAdded.x);
            viAbsDiff.y = abs(piCurrent.y-piLastAdded.y);

            if (viAbsDiff.x+viAbsDiff.y==1)
            {
                intPath.push_back(piCurrent);
                piLastAdded = piCurrent;
            }
            else {
                // Draw digital line from piLastAdded to piLastCurrent
                CIntegerPath2D pathDigitalLineSeg;
                CIntegerPath2D::const_iterator itPointSeg;

                if (viAbsDiff.x==1 && viAbsDiff.y==1)
                {
                    pathDigitalLineSeg.push_back(piLastAdded);
                    pathDigitalLineSeg.push_back(CCouple<int>(piLastAdded.x, piCurrent.y));
                    pathDigitalLineSeg.push_back(piCurrent);
                }
                else {
                    pathDigitalLineSeg.SetDigitalLineSegment4connected(piLastAdded, piCurrent);
                }

                for (itPointSeg=++pathDigitalLineSeg.begin(); itPointSeg!=pathDigitalLineSeg.end(); itPointSeg++)
                {
                    piCurrent = *itPointSeg;
                    intPath.push_back(piCurrent);
                    piLastAdded = piCurrent;
                }
            }
        }
    }

    intPath.rgbEdges = realPath.rgbEdges;
    intPath.rgbVertices = realPath.rgbVertices;
    intPath.bDisplayEdges = realPath.bDisplayVertices;
    intPath.bDisplayEdges = realPath.bDisplayEdges;
}

bool CPiecewiseGeodesicCombinationBase::MakeTempPathLists(const vector<unsigned int> &vectCombination,
    list<const CRealPath2D *> &listPtrPaths, list<const CIntegerPath2D *> &listPtrDigitalPaths) const
{
    list<list<CRealPath2D> >::const_iterator itListPaths;
    list<list<CIntegerPath2D> >::const_iterator itListDigitalPaths;
    list<CRealPath2D>::const_iterator itPath;
    list<CIntegerPath2D>::const_iterator itDigitalPath;
    const vector<unsigned int> *pCombination;
    unsigned int idxVertex, idxPath;

    listPtrPaths.clear();
    listPtrDigitalPaths.clear();

    if (vectCombination.size()==0)
        pCombination = &vectBestCombination;
    else
        pCombination = &vectCombination;

    if (pCombination->size()!=listListsAdmissiblePaths.size())
    {
        cerr<<"ERROR in CPiecewiseGeodesicCombinationBase::MakeTempPathLists():";
        cerr<<" Combination and list of admissible sets have different sizes"<<endl;
        return false;
    }

    for (idxVertex=0, itListPaths = listListsAdmissiblePaths.begin(), itListDigitalPaths = listListsAdmissiblePathsDigital.begin();
        idxVertex<listVertices.size();
        idxVertex++, itListPaths++, itListDigitalPaths++)
    {
        if (itListPaths->size()!=itListDigitalPaths->size())
        {
            cerr<<"ERROR in CPiecewiseGeodesicCombinationBase::MakeTempPathLists():";
            cerr<<" Admissible set "<<idxVertex<<": ";
            cerr<<" list of admissible paths and list of digital admissible paths have different sizes"<<endl;
            return false;
        }

        if ((*pCombination)[idxVertex]>=itListPaths->size())
        {
            cerr<<"ERROR in CPiecewiseGeodesicCombinationBase::MakeTempPathLists():";
            cerr<<" Admissible set "<<idxVertex<<": ";
            cerr<<" Label="<<vectCombination[idxVertex]<<" but set contains only "<<itListPaths->size()<<" paths"<<endl;
        }

        // Iterate through lists until the actual label is met
        for (idxPath=0, itPath=itListPaths->begin(), itDigitalPath=itListDigitalPaths->begin();
            idxPath<(*pCombination)[idxVertex];
            idxPath++, itPath++, itDigitalPath++);

        listPtrPaths.push_back(&(*itPath));
        listPtrDigitalPaths.push_back(&(*itDigitalPath));
    }

    return true;
}

bool CPiecewiseGeodesicCombinationBase::GradientDescentActionMap(
    const CCouple<float> &pfStart, const CArray2D<float> &arrayActionDerX,
    const CArray2D<float> &arrayActionDerY, CRealPath2D &path) const
{
    int iWidth, iHeight;
    CCouple<float> pfCurrent, vfGrad, pfHalf(0.5f, 0.5f);
    CCouple<int> piCurrent;
    const float *DUx, *DUy;
    int p;

	iWidth = pInputImage->GetWidth();
	iHeight = pInputImage->GetHeight();

    // Initialize a grid to make interpolation scheme available
    GeoComp::Grid8<int> grid(iHeight, iWidth);

    DUx = arrayActionDerX.GetBuffer();
    DUy = arrayActionDerY.GetBuffer();

    pfCurrent = pfStart;

    path.clear();

    piCurrent = (CCouple<int>)(pfCurrent + pfHalf);
    p = arrayActionDerX.GetOffset(piCurrent);
    int mx = grid.size()/10;

    if (DUx[p] == 0.0f && DUy[p] == 0.0f)
        return false;

    float outside = std::numeric_limits<float>::max();
    float dtn;
    float dt = 0.6f;
    path.push_back(pfCurrent);

    for (int brk = 0; brk < mx; brk++)
    {
        vfGrad.x = grid.nnInterpolation(pfCurrent.y, pfCurrent.x, DUx, outside);
        vfGrad.y = grid.nnInterpolation(pfCurrent.y, pfCurrent.x, DUy, outside);

        dtn = dt / vfGrad.L2Norm();
        pfCurrent -= vfGrad*dtn;
        //r = r - dtn * dr;
        // c = c - dtn * dc;

        if (!(pfCurrent.x>=0 && pfCurrent.x<iWidth && pfCurrent.y>=0 && pfCurrent.y<iHeight))
            return false;

        if ((path.back()-pfCurrent).L2Norm2()>=1.0f)
            path.push_back(pfCurrent);

        piCurrent = (CCouple<int>)(pfCurrent + pfHalf);
        // ri = (IT)round(r); ci = (IT)round(c);
        p = arrayActionDerX.GetOffset(piCurrent);

        if (DUx[p] == 0.0f && DUy[p] == 0.0f)
        {
            path.push_back((CCouple<float>)piCurrent);
            return true;
        }
    }

    return false;
}

float CPiecewiseGeodesicCombinationBase::GetSignedArea(const vector<unsigned int> &vectCombination) const
{
	list<const CRealPath2D *> listPtrPaths;
	list<const CIntegerPath2D *> listPtrDigitalPaths;
	list<const CIntegerPath2D *>::const_iterator itDigitalPathPtr;
	CIntegerPath2D::const_iterator itPoint;
	CCouple<int> piStart, piCurrent, piNext, piDiff;
    int iInnerRegionArea;

	MakeTempPathLists(vectCombination, listPtrPaths, listPtrDigitalPaths);

	iInnerRegionArea = 0;

    piStart = listPtrDigitalPaths.front()->front();
    piCurrent = piStart;

	for (itDigitalPathPtr=listPtrDigitalPaths.begin(); itDigitalPathPtr!=listPtrDigitalPaths.end(); itDigitalPathPtr++)
	{
	    for (itPoint=(*itDigitalPathPtr)->begin(); itPoint!=(*itDigitalPathPtr)->end(); itPoint++)
        {
            piNext = *itPoint;

            if (piNext!=piCurrent)
            {
                piDiff = piNext - piCurrent;

                if (piDiff.x==1)
                    iInnerRegionArea -= piCurrent.y;
                else if (piDiff.x==-1)
                    iInnerRegionArea += piCurrent.y;
                else if (piDiff.y==1)
                    iInnerRegionArea += piCurrent.x;
                else if (piDiff.y==-1)
                    iInnerRegionArea -= piCurrent.x;

                piCurrent = piNext;
            }
		}
	}

    if (piCurrent!=piStart)
    {
        cout<<"WARNING in CPiecewiseGeodesicCombinationBase::GetSignedArea(): "<<endl;
        cout<<"Combination of paths is not closed"<<endl;
        cout<<"  Ended at "<<piCurrent<<". Started at "<<piStart<<endl;
    }

	return (float)iInnerRegionArea*0.5f;
}

float CPiecewiseGeodesicCombinationBase::GetEnergy(const vector<unsigned int> &vectCombination)
{
	float fEnergy;

	fEnergy = 0.0f;

	if (fWeightSimplicity!=0.0f)
	{
	    GetEnergySimplicity(vectCombination);
	    if (fEnergySimplicity!=numeric_limits<float>::max())
            fEnergy += fWeightSimplicity*fEnergySimplicity;
        else
            fEnergy = numeric_limits<float>::max();
	}

    if (fWeightEdge!=0.0f)
    {
        GetEnergyEdge(vectCombination);
        if (fEnergy!=numeric_limits<float>::max())
            fEnergy += fWeightEdge*fEnergyEdge;
    }

	if (fWeightRegion!=0.0f)
	{
		GetEnergyRegion(vectCombination);
		if (fEnergyRegion!=numeric_limits<float>::max() && fEnergy!=numeric_limits<float>::max())
            fEnergy += fWeightRegion*fEnergyRegion;
        else
            fEnergy = numeric_limits<float>::max();
	}

	return fEnergy;
}

float CPiecewiseGeodesicCombinationBase::GetEnergySimplicity(const vector<unsigned int> &vectCombination)
{
    list<const CRealPath2D *> listPtrPaths;
	list<const CIntegerPath2D *> listPtrDigitalPaths;
	list<const CRealPath2D *>::const_iterator itPathPtr;
	list<const CIntegerPath2D *>::const_iterator itDigitalPathPtr, itDigitalPathPtrOther;
    CIntegerPath2D pathMerged;
    float fLength, fSelfTangency, fTwisting;

	MakeTempPathLists(vectCombination, listPtrPaths, listPtrDigitalPaths);

    fSelfTangency = 0.0f;
    fLength = 0.0f;

    // Self-tangency measure
    for (itPathPtr=listPtrPaths.begin(), itDigitalPathPtr=listPtrDigitalPaths.begin();
        itDigitalPathPtr!=listPtrDigitalPaths.end(); itPathPtr++, itDigitalPathPtr++)
    {
        itDigitalPathPtrOther = itDigitalPathPtr;
        for (++itDigitalPathPtrOther; itDigitalPathPtrOther!=listPtrDigitalPaths.end();itDigitalPathPtrOther++)
            fSelfTangency += Overlap(**itDigitalPathPtr, **itDigitalPathPtrOther);
        fLength += (*itPathPtr)->GetLength();
    }
    fEnergySimplicity = fSelfTangency/fLength;

    // Twisting measure
    pathMerged = *listPtrDigitalPaths.front();
    for (itDigitalPathPtr=++listPtrDigitalPaths.begin(); itDigitalPathPtr!=listPtrDigitalPaths.end(); itDigitalPathPtr++)
        pathMerged.AddDigitalPath4connected(**itDigitalPathPtr);
    fTwisting = InvertedAreaClosed(pathMerged);

    if (fTwisting!=numeric_limits<float>::max())
        fEnergySimplicity += fTwisting/fLength;
    else
        fEnergySimplicity = numeric_limits<float>::max();

    return fEnergySimplicity;
}

float CPiecewiseGeodesicCombinationBase::GetEnergyEdge(const vector<unsigned int> &vectCombination)
{
    list<const CRealPath2D *> listPtrPaths;
	list<const CIntegerPath2D *> listPtrDigitalPaths;
	list<const CRealPath2D *>::const_iterator itPathPtr;
	list<const CIntegerPath2D *>::const_iterator itDigitalPathPtr;
	CIntegerPath2D::const_iterator itPoint;
    float fLength, fLengthDigital, fSumGradient;

	MakeTempPathLists(vectCombination, listPtrPaths, listPtrDigitalPaths);

    fSumGradient = 0.0f;
    fLength = 0.0f;
    fLengthDigital = 0.0f;

    for (itPathPtr=listPtrPaths.begin(), itDigitalPathPtr=listPtrDigitalPaths.begin();
        itDigitalPathPtr!=listPtrDigitalPaths.end(); itPathPtr++, itDigitalPathPtr++)
    {
        for (itPoint=(*itDigitalPathPtr)->begin(); itPoint!=--(*itDigitalPathPtr)->end(); itPoint++)
        {
            fSumGradient += (1.0f - arrayGradientNorm.Element(*itPoint));
        }
        fLength += (*itPathPtr)->GetLength();
        fLengthDigital += (*itDigitalPathPtr)->GetLength();
    }
    fEnergyEdge = fSumGradient/fLengthDigital;

    return fEnergyEdge;
}

const list<CRealPath2D> &CPiecewiseGeodesicCombinationBase::GetAdmissibleSet(unsigned int idxVertex) const
{
    list<list<CRealPath2D> >::const_iterator itListPaths;
    unsigned int i;

    if (idxVertex>=listListsAdmissiblePaths.size())
    {
        cerr<<"ERROR in CPiecewiseGeodesicCombinationBase::GetAdmissibleSet(): bad index"<<endl;
        exit(-1);
    }

    itListPaths = listListsAdmissiblePaths.begin();
    for (i=0; i<idxVertex;i++)
        itListPaths++;

    return *itListPaths;
}

const CRealPath2D &CPiecewiseGeodesicCombinationBase::GetChosenPathInAdmissibleSet(unsigned int idxVertex) const
{
    list<list<CRealPath2D> >::const_iterator itListPaths;
    list<CRealPath2D>::const_iterator itPath;
    unsigned int i;

    if (idxVertex>=listListsAdmissiblePaths.size())
    {
        cerr<<"ERROR in CPiecewiseGeodesicCombinationBase::GetChosenPathInAdmissibleSet(): bad index"<<endl;
        exit(-1);
    }

    itListPaths = listListsAdmissiblePaths.begin();
    for (i=0; i<idxVertex;i++)
        itListPaths++;

    itPath = itListPaths->begin();
    for (i=0; i<vectBestCombination[idxVertex]; i++)
        itPath++;

    return *itPath;
}

CIntegerPath2D CPiecewiseGeodesicCombinationBase::GetMedialCurve(unsigned int idxVertex) const
{
    list<CVoronoiEdge>::const_iterator itVoronoiEdge;
    CVoronoiEdge::const_iterator itVoronoiPoint;
    CIntegerPath2D pathMedial;
    unsigned int i;

    if (idxVertex>=listVoronoiEdges.size())
    {
        cerr<<"ERROR in CPiecewiseGeodesicCombinationBase::GetMedialCurve(): bad index"<<endl;
        exit(-1);
    }

    itVoronoiEdge = listVoronoiEdges.begin();
    for (i=0; i<idxVertex; i++)
        itVoronoiEdge++;

    pathMedial = itVoronoiEdge->ToIntegerPath();
    pathMedial.rgbEdges.Set(0, 0, 0);

    return pathMedial;
}

/**********************************************************************
*               CPiecewiseGeodesicCombinationGrayscale                *
**********************************************************************/

CPiecewiseGeodesicCombinationGrayscale::CPiecewiseGeodesicCombinationGrayscale():CPiecewiseGeodesicCombinationBase()
{
    iNbBinsPerComponent = 256;
}

void CPiecewiseGeodesicCombinationGrayscale::AttachImage(const CImage2D *pImage)
{
    CArray2D<float> arrayImage;
	int iBitsPerPixel;
	CArray2D<float> arrayGaussian, arrayDiffX, arrayDiffY;
    CImage2D imageGradient;
    const float *pDiffX, *pDiffY;
    float *pNorm;
    const int *pImageScaled;
    int i, iSize;

	iBitsPerPixel = pImage->GetBitsPerPixel();
	if (iBitsPerPixel!=8)
	{
		cerr<<"ERROR in CPiecewiseGeodesicCombinationGrayscale::AttachImage(): instances of CPiecewiseGeodesicCombinationGrayscale can be attached to 8-bit images only"<<endl;
		pInputImage = NULL;
		return;
	}

	pInputImage = pImage;
	pInputImage->ConvertToArray2DFloat(arrayImage);

    histoInside.SetRange(0.0f, 1.0f);
    histoOutside.SetRange(0.0f, 1.0f);
    histoImage.SetRange(0.0f, 1.0f);

    histoInside.SetNbBins(iNbBinsPerComponent);
    histoOutside.SetNbBins(iNbBinsPerComponent);
    histoImage.SetNbBins(iNbBinsPerComponent);

    histoInside.SetGaussianStdDeviation(2.0f);
    histoOutside.SetGaussianStdDeviation(2.0f);
    histoImage.SetGaussianStdDeviation(2.0f);

    // Add offset to the histogram scaled values
    histoInside.MakeImageScaled(arrayImage, arrayImageHistogramScaled);

    // Compute histogram of entire image
    pImageScaled = arrayImageHistogramScaled.GetBuffer();
    iSize = arrayImageHistogramScaled.GetWidth()*arrayImageHistogramScaled.GetHeight();
    for (i=0; i<iSize; i++)
    {
        histoImage.AddElementScaled(*pImageScaled);
        pImageScaled++;
    }

    // Compute gradient magnitude
    // Gaussian mask's half size should be at least 3
    if (fGradientScale<1.0f)
        arrayGaussian.SetGaussianKernel(fGradientScale, 3);
    else
        arrayGaussian.SetGaussianKernel(fGradientScale);

    // Convolve the image with x and y-derivatives of gaussian
    arrayDiffX = arrayImage.Convolve(arrayGaussian.DerivativeX(1, ARRAYNDFLOAT_CENTERED));
    arrayDiffY = arrayImage.Convolve(arrayGaussian.DerivativeY(1, ARRAYNDFLOAT_CENTERED));

    // Compute the magnitude of the obtained "smooth gradient" for every pixel
    arrayGradientNorm.Init(arrayImage.GetSize());
    pDiffX = arrayDiffX.GetBuffer();
    pDiffY = arrayDiffY.GetBuffer();
    pNorm = arrayGradientNorm.GetBuffer();
    iSize = arrayGradientNorm.GetWidth()*arrayGradientNorm.GetHeight();
    for (i=0; i<iSize; i++)
    {
        *pNorm = sqrt((*pDiffX)*(*pDiffX) + (*pDiffY)*(*pDiffY));
        pDiffX++;
        pDiffY++;
        pNorm++;
    }

	UpdatePotential();
}

void CPiecewiseGeodesicCombinationGrayscale::UpdateHistograms(const vector<unsigned int> &vectCombination)
{
    list<const CRealPath2D *> listPtrPaths;
	list<const CIntegerPath2D *> listPtrDigitalPaths;
	list<const CIntegerPath2D *>::const_iterator itDigitalPathPtr;
	CIntegerPath2D::const_iterator itPoint;
	CCouple<int> piIntegral, piStart, piCurrent, piNext, piDiff;

	MakeTempPathLists(vectCombination, listPtrPaths, listPtrDigitalPaths);

	histoInside.Empty();
	histoOutside.Empty();

    piStart = listPtrDigitalPaths.front()->front();
    piCurrent = piStart;

	for (itDigitalPathPtr=listPtrDigitalPaths.begin(); itDigitalPathPtr!=listPtrDigitalPaths.end(); itDigitalPathPtr++)
	{
	    for (itPoint=(*itDigitalPathPtr)->begin(); itPoint!=(*itDigitalPathPtr)->end(); itPoint++)
        {
            piNext = *itPoint;

            if (piNext!=piCurrent)
            {
                piDiff = piNext - piCurrent;

                if (piDiff.x==1)
                {
                    piIntegral.x = piCurrent.x;
                    for (piIntegral.y=0;piIntegral.y<piCurrent.y;piIntegral.y++)
                        histoInside.SubtractElementScaled(arrayImageHistogramScaled.Element(piIntegral), 0.5f);
                }
                else if (piDiff.x==-1)
                {
                    piIntegral.x = piCurrent.x-1;
                    for (piIntegral.y=0;piIntegral.y<piCurrent.y;piIntegral.y++)
                        histoInside.AddElementScaled(arrayImageHistogramScaled.Element(piIntegral), 0.5f);
                }
                else if (piDiff.y==1)
                {
                    piIntegral.y = piCurrent.y;
                    for (piIntegral.x=0;piIntegral.x<piCurrent.x;piIntegral.x++)
                        histoInside.AddElementScaled(arrayImageHistogramScaled.Element(piIntegral), 0.5f);
                }
                else if (piDiff.y==-1)
                {
                    piIntegral.y = piCurrent.y-1;
                    for (piIntegral.x=0;piIntegral.x<piCurrent.x;piIntegral.x++)
                        histoInside.SubtractElementScaled(arrayImageHistogramScaled.Element(piIntegral), 0.5f);
                }

                piCurrent = piNext;
            }
		}
	}

    if (piCurrent!=piStart)
    {
        cout<<"WARNING in CPiecewiseGeodesicCombinationGrayscale::UpdateHistograms(...): Combination of paths is not closed"<<endl;
        cout<<"  Ended at "<<piCurrent<<". Started at "<<piStart<<endl;
    }

    histoOutside = histoImage - histoInside;

    /*
    CArray2D<float> arrayHistoTemp(histoInside.GetSize(), 1);

    #ifdef GRAPHICWND_H
    static CArray2DFloatWnd *pWndHistoInside = NULL;
    if (pWndHistoInside==NULL)
        pWndHistoInside = new CArray2DFloatWnd(256, 50, "Histo inside");

    for (int x=0; x<histoInside.GetSize(); x++)
        arrayHistoTemp.Element(x,0) = histoOutside[x];

    pWndHistoInside->SetArray(arrayHistoTemp);
    pWndHistoInside->UpdateImage();
    pWndHistoInside->show();
    #endif
    */
    histoInside.GaussianSmooth();
	histoOutside.GaussianSmooth();
    /*
    #ifdef GRAPHICWND_H
    static CArray2DFloatWnd *pWndHistoInsideSmoothed = NULL;
    if (pWndHistoInsideSmoothed==NULL)
        pWndHistoInsideSmoothed = new CArray2DFloatWnd(256, 50, "Histo inside smoothed");

    for (int x=0; x<histoInside.GetSize(); x++)
        arrayHistoTemp.Element(x,0) = histoOutside[x];

    pWndHistoInsideSmoothed->SetArray(arrayHistoTemp);
    pWndHistoInsideSmoothed->UpdateImage();
    pWndHistoInsideSmoothed->show();
    #endif
    */
}

float CPiecewiseGeodesicCombinationGrayscale::GetEnergyRegion(const vector<unsigned int> &vectCombination)
{
    UpdateHistograms(vectCombination);

    if (histoInside.GetSumWeights()<0.0f)
        fEnergyRegion = numeric_limits<float>::max();
    else
        fEnergyRegion = histoInside.BhattacharyyaCoefficient(histoOutside);

    return fEnergyRegion;
}

/**********************************************************************
*                 CPiecewiseGeodesicCombinationColor                  *
**********************************************************************/

CPiecewiseGeodesicCombinationColor::CPiecewiseGeodesicCombinationColor():CPiecewiseGeodesicCombinationBase()
{
    // Default color space is RGB
	iColorSpace = COLORSPACE_RGB;

	bIgnoreBrightnessComponent = false;

	iNbBinsPerComponent = 32;
}

void CPiecewiseGeodesicCombinationColor::AttachImage(const CImage2D *pImage)
{
    CArray2D<CTriplet<float> > arrayImage;
	int iBitsPerPixel;
	CTriplet<float> fHistoRangeMin(0.0f, 0.0f, 0.0f), fHistoRangeMax(1.0f, 1.0f, 1.0f);
    CTriplet<int> iNbBins;
    CArray2D<CTriplet<float> > arrayDiffX, arrayDiffY;
    CArray2D<float> arrayGaussian;
    // CImage2D imageGradient;
    const CTriplet<float> *pDiffX, *pDiffY;
    const CTriplet<int> *pImageScaled;
    float *pNorm;
    int i, iSize;

	iBitsPerPixel = pImage->GetBitsPerPixel();
	if (iBitsPerPixel!=24 && iBitsPerPixel!=32)
	{
		cerr<<"ERROR in CPiecewiseGeodesicCombinationColor::AttachImage(): instances of CPiecewiseGeodesicCombinationColor can be attached to RGB images only"<<endl;
		pInputImage = NULL;
		return;
	}

	pInputImage = pImage;

	if (iColorSpace==COLORSPACE_RGB)
	{
		pInputImage->ConvertToArray2DTripletFloatRGB(arrayImage);
        fHistoRangeMin.Set(0.0f, 0.0f, 0.0f);
        fHistoRangeMax.Set(1.0f, 1.0f, 1.0f);
	}
	else if (iColorSpace==COLORSPACE_YUV)
	{
        pInputImage->ConvertToArray2DTripletFloatYUV(arrayImage);
        fHistoRangeMin.Set(0.0f, -0.5f, -0.65f);
        fHistoRangeMax.Set(1.0f, 0.5f, 0.65f);
	}
	else if (iColorSpace==COLORSPACE_LAB)
	{
		pInputImage->ConvertToArray2DTripletFloatLAB(arrayImage);
        fHistoRangeMin.Set(0.0f, -0.55f, -0.5f);
        fHistoRangeMax.Set(1.0f, 0.7f, 0.55f);
    }

	// If chosen color space separates brightness from color information (i.e. not RGB)
	// and brightness component should be ignored for illumination invariance purpose,
	// then set brightness component to zero everywhere
	if (iColorSpace!=COLORSPACE_RGB && bIgnoreBrightnessComponent==true)
	{
		CTriplet<float> *pColor;
		int i,iSize;

		iSize = arrayImage.GetWidth()*arrayImage.GetHeight();
		pColor = arrayImage.GetBuffer();
		for (i=0; i<iSize; i++)
		{
			pColor->x = 0.0f;
			pColor++;
		}

        iNbBins.Set(1,iNbBinsPerComponent,iNbBinsPerComponent);
	}
	else
        iNbBins.Set(iNbBinsPerComponent, iNbBinsPerComponent, iNbBinsPerComponent);

    histoInside.SetRange(fHistoRangeMin, fHistoRangeMax);
    histoOutside.SetRange(fHistoRangeMin, fHistoRangeMax);
    histoImage.SetRange(fHistoRangeMin, fHistoRangeMax);

    histoInside.SetNbBins(iNbBins);
    histoOutside.SetNbBins(iNbBins);
    histoImage.SetNbBins(iNbBins);

	histoInside.SetGaussianStdDeviation(2.0f);
    histoOutside.SetGaussianStdDeviation(2.0f);
    histoImage.SetGaussianStdDeviation(2.0f);

    // Add offset to the histogram scaled values
    histoInside.MakeImageScaled(arrayImage, arrayImageHistogramScaled);

    // Compute histogram of entire image
    pImageScaled = arrayImageHistogramScaled.GetBuffer();
    iSize = arrayImageHistogramScaled.GetWidth()*arrayImageHistogramScaled.GetHeight();
    for (i=0; i<iSize; i++)
    {
        histoImage.AddElementScaled(*pImageScaled);
        pImageScaled++;
    }

	// Compute gradient magnitude
    // Gaussian mask's half size should be at least 3
    if (fGradientScale<1.0f)
        arrayGaussian.SetGaussianKernel(fGradientScale, 3);
    else
        arrayGaussian.SetGaussianKernel(fGradientScale);

    // Convolve the image with x and y-derivatives of gaussian
    arrayDiffX = arrayImage.Convolve(arrayGaussian.DerivativeX(1, ARRAYNDFLOAT_CENTERED));
    arrayDiffY = arrayImage.Convolve(arrayGaussian.DerivativeY(1, ARRAYNDFLOAT_CENTERED));

    // Compute the magnitude of the obtained "smooth gradient" for every pixel
    arrayGradientNorm.Init(arrayImage.GetSize());
    pDiffX = arrayDiffX.GetBuffer();
    pDiffY = arrayDiffY.GetBuffer();
    pNorm = arrayGradientNorm.GetBuffer();
    iSize = arrayGradientNorm.GetWidth()*arrayGradientNorm.GetHeight();
    for (i=0; i<iSize; i++)
    {
        *pNorm = sqrt(pDiffX->L2Norm2() + pDiffY->L2Norm2());
        pDiffX++;
        pDiffY++;
        pNorm++;
    }

	UpdatePotential();
}

void CPiecewiseGeodesicCombinationColor::UpdateHistograms(const vector<unsigned int> &vectCombination)
{
    list<const CRealPath2D *> listPtrPaths;
	list<const CIntegerPath2D *> listPtrDigitalPaths;
	list<const CIntegerPath2D *>::const_iterator itDigitalPathPtr;
	CIntegerPath2D::const_iterator itPoint;
	CCouple<int> piIntegral, piStart, piCurrent, piNext, piDiff;

	MakeTempPathLists(vectCombination, listPtrPaths, listPtrDigitalPaths);

	histoInside.Empty();
	histoOutside.Empty();

    piStart = listPtrDigitalPaths.front()->front();
    piCurrent = piStart;

	for (itDigitalPathPtr=listPtrDigitalPaths.begin(); itDigitalPathPtr!=listPtrDigitalPaths.end(); itDigitalPathPtr++)
	{
	    for (itPoint=(*itDigitalPathPtr)->begin(); itPoint!=(*itDigitalPathPtr)->end(); itPoint++)
        {
            piNext = *itPoint;

            if (piNext!=piCurrent)
            {
                piDiff = piNext - piCurrent;

                if (piDiff.x==1)
                {
                    piIntegral.x = piCurrent.x;
                    for (piIntegral.y=0;piIntegral.y<piCurrent.y;piIntegral.y++)
                        histoInside.SubtractElementScaled(arrayImageHistogramScaled.Element(piIntegral), 0.5f);
                    //fInnerRegionArea -= (float)piCurrent.y;
                }
                else if (piDiff.x==-1)
                {
                    piIntegral.x = piCurrent.x-1;
                    for (piIntegral.y=0;piIntegral.y<piCurrent.y;piIntegral.y++)
                        histoInside.AddElementScaled(arrayImageHistogramScaled.Element(piIntegral), 0.5f);
                    // fInnerRegionArea += (float)piCurrent.y;
                }
                else if (piDiff.y==1)
                {
                    piIntegral.y = piCurrent.y;
                    for (piIntegral.x=0;piIntegral.x<piCurrent.x;piIntegral.x++)
                        histoInside.AddElementScaled(arrayImageHistogramScaled.Element(piIntegral), 0.5f);
                    //fInnerRegionArea += (float)piCurrent.x;
                }
                else if (piDiff.y==-1)
                {
                    piIntegral.y = piCurrent.y-1;
                    for (piIntegral.x=0;piIntegral.x<piCurrent.x;piIntegral.x++)
                        histoInside.SubtractElementScaled(arrayImageHistogramScaled.Element(piIntegral), 0.5f);
                    // fInnerRegionArea -= (float)piCurrent.x;
                }

                piCurrent = piNext;
            }
		}
	}

    if (piCurrent!=piStart)
    {
        cout<<"WARNING in CPiecewiseGeodesicCombinationColor::UpdateHistograms(...): Combination of paths is not closed"<<endl;
        cout<<"  Ended at "<<piCurrent<<". Started at "<<piStart<<endl;
    }

	// fInnerRegionArea*=0.5f;
	// fOuterRegionArea = (float)(iWidth*iHeight)-fInnerRegionArea;

    histoOutside = histoImage - histoInside;

    histoInside.GaussianSmooth();
	histoOutside.GaussianSmooth();

	#ifdef GRAPHICWND_H
	char strWindowName[100], strNumber[10], strCombination[50]="";
	unsigned int idxVertex;

	if (vectCombination.size()==0)
    {
        for (idxVertex=0; idxVertex<vectBestCombination.size(); idxVertex++)
        {
            sprintf(strNumber, "%d ", vectBestCombination[idxVertex]);
            strcat(strCombination, strNumber);
        }
    }
    else {
        for (idxVertex=0; idxVertex<vectCombination.size(); idxVertex++)
        {
            sprintf(strNumber, "%d ", vectCombination[idxVertex]);
            strcat(strCombination, strNumber);
        }
    }

    sprintf(strWindowName, "Histo inside %s", strCombination);
    CArray3DFloatWnd::Show(histoInside, strWindowName);

    sprintf(strWindowName, "Histo outside %s", strCombination);
    CArray3DFloatWnd::Show(histoOutside, strWindowName);

    CArray3DFloatWnd::Show(histoImage, "Histo image");
	#endif
}

float CPiecewiseGeodesicCombinationColor::GetEnergyRegion(const vector<unsigned int> &vectCombination)
{
    UpdateHistograms(vectCombination);

    if (histoInside.GetSumWeights()<0.0f || histoOutside.GetSumWeights()<0.0f)
        fEnergyRegion = numeric_limits<float>::max();
    else
        fEnergyRegion = histoInside.BhattacharyyaCoefficient(histoOutside);

    return fEnergyRegion;
}

