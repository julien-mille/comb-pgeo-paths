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
path.cpp

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

#include <vector>
#include "path.h"


/*********************************************************************
*                        CIntegerPath2D                              *
*********************************************************************/

CIntegerPath2D::CIntegerPath2D(const CRealPath2D &realPath)
{
	list<CCouple<float> >::const_iterator itPoint;
	CCouple<float> pfHalf(0.5f, 0.5f);

	clear();
	for (itPoint=realPath.begin(); itPoint!=realPath.end(); itPoint++)
		push_back((CCouple<int>)(*itPoint + pfHalf));
}

float CIntegerPath2D::GetLength() const
{
	list<CCouple<int> >::const_iterator itPoint, itPointNext;
	float fLength;

	fLength = 0.0f;
	for (itPoint=begin(); itPoint!=--end(); itPoint++)
	{
	    itPointNext = itPoint;
	    itPointNext++;

		fLength += sqrt((float)((*itPointNext - *itPoint).L2Norm2()));
	}
	return fLength;
}

void CIntegerPath2D::SetDigitalLineSegment4connected(const CCouple<int> &piStart, const CCouple<int> &piEnd)
{
    CCouple<int> piCurrent;
    float fX, fY, a;
	int iStep, iY, iX;

	clear();

	if (piStart==piEnd)
		push_back(piStart);
	else if (piStart.x!=piEnd.x)
	{
		a = (float)(piEnd.y-piStart.y)/(piEnd.x-piStart.x);

		if (fabs(a)<=1.0f)
		{
			// Slope is greater in x than in y
			// Loop on x and compute corresponding y
			if (piEnd.x>piStart.x)
                iStep = 1;
            else {
                iStep = -1;
                a*=-1.0f;
            }

            piCurrent = piStart;
            fY = (float)piStart.y;
			for (piCurrent.x=piStart.x; piCurrent.x!=piEnd.x+iStep; piCurrent.x+=iStep)
			{
			    push_back(piCurrent);
			    fY += a;
			    iY = (int)(fY+0.5f);
			    if (iY!=piCurrent.y)
                    push_back(CCouple<int>(piCurrent.x, iY));
                piCurrent.y = iY;
			}
		}
		else {
			// Slope is greater in y than in x
			// Loop on y and compute corresponding x
            a = 1.0f/a;

            if (piEnd.y>piStart.y)
                iStep = 1;
            else {
                iStep = -1;
                a*=-1.0f;
            }

            piCurrent = piStart;
            fX = (float)piStart.x;
			for (piCurrent.y=piStart.y; piCurrent.y!=piEnd.y+iStep; piCurrent.y+=iStep)
			{
			    push_back(piCurrent);
			    fX += a;
			    iX = (int)(fX+0.5f);
			    if (iX!=piCurrent.x)
                    push_back(CCouple<int>(iX, piCurrent.y));
			    piCurrent.x = iX;
			}
		}
	}
	else {
		if (piEnd.y>piStart.y)
            iStep = 1;
        else
            iStep = -1;

		piCurrent = piStart;
        for (piCurrent.y=piStart.y; piCurrent.y!=piEnd.y+iStep; piCurrent.y+=iStep)
		    push_back(piCurrent);
	}
}

void CIntegerPath2D::AddDigitalPath4connected(const CIntegerPath2D &path)
{
    CIntegerPath2D::const_iterator itPoint, itBegin, itEnd;
    CCouple<int> viDiff;

    // If one of the two paths is empty
    if (size()==0)
    {
        *this = path;
        return;
    }
    else if (path.size()==0)
        return;

    viDiff = path.front() - back();

    if (viDiff.x==0 && viDiff.y==0)
        itPoint = ++path.begin();
    else if (abs(viDiff.x)+abs(viDiff.y)==1)
        itPoint = path.begin();
    else {
        cerr<<"ERROR in CIntegerPath2D::AddDigitalPath4connected(): end-points do not match"<<endl;
        cerr<<"Last point of current path = "<<back()<<". First point of path to add = "<<path.front()<<endl;
        return;
    }

    for (; itPoint!=path.end(); itPoint++)
        push_back(*itPoint);

    if (back()==front())
        pop_back();
}

float CIntegerPath2D::GetArea4connected() const
{
    CIntegerPath2D::const_iterator itPoint;
    CCouple<int> piCurrent, piNext, piDiff;
    int iArea;

    piCurrent = front();
    iArea = 0;
    for (itPoint=++begin(); itPoint!=end(); itPoint++)
    {
        piNext = *itPoint;

        piDiff = piNext - piCurrent;
        if (abs(piDiff.x) + abs(piDiff.y)<=1)
        {
            if (piDiff.x==1)
                iArea -= piCurrent.y;
            else if (piDiff.x==-1)
                iArea += piCurrent.y;
            else if (piDiff.y==1)
                iArea += piCurrent.x;
            else if (piDiff.y==-1)
                iArea -= piCurrent.x;
            piCurrent = piNext;
        }
        else {
            cerr<<"ERROR in CIntegerPath2D::GetArea4connected(): ";
            cerr<<" path is not 4-connected between "<<piCurrent<<" and "<<piNext<<endl;
            return 0.0f;
        }
    }

    piCurrent = back();
    piNext = front();

    if (piNext!=piCurrent)
    {
        piDiff = piNext - piCurrent;

        if (abs(piDiff.x) + abs(piDiff.y)<=1)
        {
            if (piDiff.x==1)
                iArea -= piCurrent.y;
            else if (piDiff.x==-1)
                iArea += piCurrent.y;
            else if (piDiff.y==1)
                iArea += piCurrent.x;
            else if (piDiff.y==-1)
                iArea -= piCurrent.x;
        }
        else {
            cerr<<"ERROR in CIntegerPath2D::GetArea4connected(): ";
            cerr<<" path is not closed between "<<piCurrent<<" and "<<piNext<<endl;
            return 0.0f;
        }
    }

    return (float)iArea*0.5f;
}

void CIntegerPath2D::DrawInImage(CImage2D &image2D, int iZoom) const
{
    list<CCouple<int> >::const_iterator itPoint;
    CArray1D<CCouple<int> > vectIntPoints;
    unsigned int i;

    if (size()==0)
        return;

    vectIntPoints.Init(size());

    for (itPoint=begin(), i=0; itPoint!=end(); itPoint++, i++)
        vectIntPoints[i] = (*itPoint)*iZoom + CCouple<int>(iZoom, iZoom)/2; // - CCouple<int>(1,1);

	if (image2D.GetBitsPerPixel()==8)
	{
		CImage2DPixel col;

		if (bDisplayEdges==true)
		{
		    col = rgbEdges.ToGrayScale();
		    for (i=0;i<size()-1;i++)
				image2D.DrawLine(vectIntPoints[i], vectIntPoints[i+1], col);
		}
		/*
		if (bDisplayVertices==true)
		{
			col = rgbVertices.ToGrayScale();
			for (i=0; i<size(); i++)
				image2D.DrawFilledRectangle(vectIntPoints[i]-CCouple<int>(1,1), vectIntPoints[i]+CCouple<int>(1,1), col);
		}*/
	}
	else if (image2D.GetBitsPerPixel()==24 || image2D.GetBitsPerPixel()==32)
	{
		if (bDisplayEdges==true)
		{
		    CImage2DFloatRGBPixel fRgbStart, fRgbEnd, fRgbMix;
		    float fCoef;

		    fRgbStart = (CImage2DFloatRGBPixel)rgbEdges;
		    fRgbEnd = CImage2DFloatRGBPixel(255.0f, 255.0f, 255.0f); // - fRgbStart;

			for (i=0;i<size()-1;i++)
			{
			    fCoef = (float)i/(float)size();
			    fRgbMix = fRgbStart*(1.0f-fCoef) + fRgbEnd*fCoef;
			    image2D.DrawLineRGB(vectIntPoints[i], vectIntPoints[i+1], fRgbMix.ToByteRGBPixel());
			}
		}
		/*
		if (bDisplayVertices==true)
		{
			for (i=0; i<size(); i++)
				image2D.DrawFilledRectangleRGB(vectIntPoints[i]-CCouple<int>(1,1), vectIntPoints[i]+CCouple<int>(1,1), rgbVertices);
		}*/
	}
}


void CIntegerPath2D::DrawInImagePixelCorners(CImage2D &image2D, int iZoom) const
{
    list<CCouple<int> >::const_iterator itPoint;
    CArray1D<CCouple<int> > vectIntPoints;
    unsigned int i;

    if (size()==0)
        return;

    vectIntPoints.Init(size());

    for (itPoint=begin(), i=0; itPoint!=end(); itPoint++, i++)
        vectIntPoints[i] = (*itPoint)*iZoom; // - CCouple<int>(1,1);

	if (image2D.GetBitsPerPixel()==8)
	{
		// int iVertex;
		// CCouple<int> *pVertex;
		// CTriplet<unsigned char> rgbVertices, rgbEdges;
		CImage2DPixel col;

		if (bDisplayEdges==true)
		{
		    col = rgbEdges.ToGrayScale();
		    for (i=0;i<size()-1;i++)
				image2D.DrawLine(vectIntPoints[i], vectIntPoints[i+1], col);
		}
		/*
		if (bDisplayVertices==true)
		{
			col = rgbVertices.ToGrayScale();
			for (i=0; i<size(); i++)
				image2D.DrawFilledRectangle(vectIntPoints[i]-CCouple<int>(1,1), vectIntPoints[i]+CCouple<int>(1,1), col);
		}*/
	}
	else if (image2D.GetBitsPerPixel()==24 || image2D.GetBitsPerPixel()==32)
	{
		if (bDisplayEdges==true)
		{
		    CImage2DFloatRGBPixel fRgbStart, fRgbEnd, fRgbMix;
		    float fCoef;

		    fRgbStart = (CImage2DFloatRGBPixel)rgbEdges;
		    fRgbEnd = CImage2DFloatRGBPixel(255.0f, 255.0f, 255.0f) - fRgbStart;

			for (i=0;i<size()-1;i++)
			{
			    fCoef = (float)i/(float)size();
			    fRgbMix = fRgbStart*(1.0f-fCoef) + fRgbEnd*fCoef;
			    image2D.DrawLineRGB(vectIntPoints[i], vectIntPoints[i+1], fRgbMix.ToByteRGBPixel());
			}
		}
		/*
		if (bDisplayVertices==true)
		{
			for (i=0; i<size(); i++)
				image2D.DrawFilledRectangleRGB(vectIntPoints[i]-CCouple<int>(1,1), vectIntPoints[i]+CCouple<int>(1,1), rgbVertices);
		}*/
	}
}

/*********************************************************************
*                           CRealPath2D                              *
*********************************************************************/

CRealPath2D::CRealPath2D(const CIntegerPath2D &intPath)
{
	list<CCouple<int> >::const_iterator itPoint;

	clear();
	for (itPoint=intPath.begin(); itPoint!=intPath.end(); itPoint++)
		push_back((CCouple<float>)(*itPoint));
}

float CRealPath2D::GetLength() const
{
	list<CCouple<float> >::const_iterator itPoint, itPointNext;
	float fLength;

	fLength = 0.0f;
	for (itPoint=begin(); itPoint!=--end(); itPoint++)
	{
	    itPointNext = itPoint;
	    itPointNext++;

		fLength += (*itPointNext - *itPoint).L2Norm();
	}
	return fLength;
}

void CRealPath2D::InitLine(const CCouple<float> &pfStart, const CCouple<float> &pfEnd)
{
	float fDist; //, fCoef, fCoefStep;
	CCouple<float> vfDir, pfPathPoint;
	int iPoint, iNbPoints;
	clear();

	fDist = (pfEnd-pfStart).L2Norm();
	vfDir = (pfEnd-pfStart)/fDist;
	// fCoefStep = 1.0f/fDist;
	iNbPoints = (int)fDist;
	pfPathPoint = pfStart;

	for (iPoint=0; iPoint<iNbPoints; iPoint++)
	{
		push_back(pfPathPoint);
		pfPathPoint += vfDir;
	}

	push_back(pfEnd);
}

void CRealPath2D::LaplacianSmooth()
{
	CArray1D<CCouple<float> > vectNewVertices;
	list<CCouple<float> >::iterator itPointPrev, itPoint, itPointNext;
	int i;

	if (size()<3)
		return;

	vectNewVertices.Init(size()-2);

	for (itPoint=++begin(), i=0; itPoint!=--end(); itPoint++, i++)
	{
	    itPointPrev = itPoint;
	    itPointPrev--;
	    itPointNext = itPoint;
	    itPointNext++;

		vectNewVertices[i] = (*itPoint + (*itPointPrev + *itPointNext)*0.5f)*0.5f;
	}

	for (itPoint=++begin(), i=0; itPoint!=--end(); itPoint++, i++)
		*itPoint = vectNewVertices[i];
}

void CRealPath2D::Translate(const CCouple<float> &vfTrans)
{
    list<CCouple<float> >::iterator itPoint;

    for (itPoint=begin(); itPoint!=end(); itPoint++)
	    *itPoint += vfTrans;
}

float CRealPath2D::GetSignedAreaWithLineSegment() const
{
    list<CCouple<float> >::const_iterator itPoint, itPointNext;
    float fArea;

    if (size()<2)
        cout<<"WARNING in CRealPath2D::GetSignedAreaWithLineSegment()"<<endl;

    fArea = 0.0f;

    for (itPoint=begin(), itPointNext=++begin(); itPointNext!=end(); itPoint++, itPointNext++)
        fArea += itPoint->x*itPointNext->y - itPointNext->x*itPoint->y;

    itPoint = --end();
    itPointNext = begin();
    fArea += itPoint->x*itPointNext->y - itPointNext->x*itPoint->y;

    return fArea*0.5f;
}

CRealPath2D CRealPath2D::Resampled(float fDist) const
{
    list<CCouple<float> >::const_iterator itPoint, itPointNext;
    CCouple<float> pfLastPointAdded;
    CRealPath2D path;
    float fDistTotal;
    float fDist2 = fDist*fDist;

    path.push_back(front());
    pfLastPointAdded = front();

    for (itPoint=begin(), itPointNext=++begin(); itPointNext!=end(); itPoint++, itPointNext++)
    {
        fDistTotal = (*itPointNext-pfLastPointAdded).L2Norm2();
        if (fDistTotal>=fDist2)
        {
            path.push_back(*itPointNext);
            pfLastPointAdded = *itPointNext;
        }
    }
    if (*itPoint!=path.back())
        path.push_back(*itPoint);

    path.rgbEdges = rgbEdges;
    path.rgbVertices = rgbVertices;
    path.bDisplayEdges = bDisplayEdges;
    path.bDisplayVertices = bDisplayVertices;

    return path;
}

void CRealPath2D::DrawInImage(CImage2D &image2D, int iZoom) const
{
    list<CCouple<float> >::const_iterator itPoint;
    CArray1D<CCouple<int> > vectIntPoints;
    unsigned int i;

    if (size()==0)
        return;

    vectIntPoints.Init(size());

    for (itPoint=begin(), i=0; itPoint!=end(); itPoint++, i++)
        vectIntPoints[i] = (CCouple<int>)((*itPoint)*(float)iZoom) + CCouple<int>(iZoom/2, iZoom/2);

	if (image2D.GetBitsPerPixel()==8)
	{
		CImage2DPixel col;

		if (bDisplayEdges==true)
		{
		    col = rgbEdges.ToGrayScale();
		    for (i=0;i<size()-1;i++)
				image2D.DrawLine(vectIntPoints[i], vectIntPoints[i+1], col);
		}
		if (bDisplayVertices==true)
		{
			col = rgbVertices.ToGrayScale();
			for (i=0; i<size(); i++)
				image2D.DrawFilledRectangle(vectIntPoints[i]-CCouple<int>(1,1), vectIntPoints[i]+CCouple<int>(1,1), col);
		}
	}
	else if (image2D.GetBitsPerPixel()==24 || image2D.GetBitsPerPixel()==32)
	{
		if (bDisplayEdges==true)
		{
			for (i=0;i<size()-1;i++)
				image2D.DrawLineRGB(vectIntPoints[i], vectIntPoints[i+1], rgbEdges);
		}
		if (bDisplayVertices==true)
		{
			for (i=0; i<size(); i++)
				image2D.DrawFilledRectangleRGB(vectIntPoints[i]-CCouple<int>(1,1), vectIntPoints[i]+CCouple<int>(1,1), rgbVertices);
		}
	}
}
