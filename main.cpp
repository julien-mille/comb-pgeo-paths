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

#include <iostream>
#include <stdio.h>

#include "piecewisegeodesiccombination.h"

using namespace std;

int main(int argc, char *argv[])
{
	CImage2D imgInput, imgOutput, imgContour;
	CPiecewiseGeodesicCombinationBase *pPGComb;
	int iZoom;

    
	CImage2D::StartImage2D();

	// Load the image file
	// if (imgInput.Load("./hand.jpg")==false)
	if (imgInput.Load("./leaf.jpg")==false)
	{
		cerr<<"ERROR: unable to open input image."<<endl;
		return -1;
	}

	// Set zooming factor for output images
	iZoom = 1;

	// Allocate the model to a grayscale or color instance depending on image format
	if (imgInput.GetBitsPerPixel()==8)
	{
        pPGComb = new CPiecewiseGeodesicCombinationGrayscale();

		// Convert input image to 24-bit RGB format for output
		imgInput.Convert8bitsto24bits(imgOutput);

		// Rescale output image if needed
		if (iZoom>1)
			imgOutput = imgOutput.ResizeRGB(imgOutput.GetWidth()*iZoom, imgOutput.GetHeight()*iZoom);
	}
	else if (imgInput.GetBitsPerPixel()==24 || imgInput.GetBitsPerPixel()==32)
	{
        CPiecewiseGeodesicCombinationColor *pPGCombColor;

        pPGComb = new CPiecewiseGeodesicCombinationColor();

        pPGCombColor = (CPiecewiseGeodesicCombinationColor *)pPGComb;

        // Set parameters specific to the color model
        pPGCombColor->iColorSpace = CPiecewiseGeodesicCombinationColor::COLORSPACE_RGB; // or COLORSPACE_RGB, COLORSPACE_LAB
        pPGCombColor->bIgnoreBrightnessComponent = false;

		// Copy the input image to the output image (with rescaling if needed)
		if (iZoom>1)
			imgOutput = imgInput.ResizeRGB(imgInput.GetWidth()*iZoom, imgInput.GetHeight()*iZoom);
		else
			imgOutput = imgInput;
	}
	else {
		cerr<<"ERROR: number of bits/pixel is not supported."<<endl;
		return -1;
	}

	// Other parameters may be modified here, before attaching the model to the image
	// ...

	// Attach model to image data
	pPGComb->AttachImage(&imgInput);

	// Initial vertices and gradient weight for 'hand.jpg'
	// pPGComb->fWeightGradientPotential = 15.0f;
	// pPGComb->AddVertex(CCouple<float>(392.0, 211.0));
	// pPGComb->AddVertex(CCouple<float>(436.0, 433.0));
    
	// Initial vertices and gradient weight for 'leaf.jpg'
	pPGComb->fWeightGradientPotential = 10.0f;
	pPGComb->AddVertex(CCouple<float>(295.0, 180.0));
	pPGComb->AddVertex(CCouple<float>(439.0, 324.0));
	pPGComb->AddVertex(CCouple<float>(253.0, 409.0));
	
	// Update potential (gradient weight was modified)
	pPGComb->UpdatePotential();

        // Build piecewise-geodesic contour
	pPGComb->UpdateAllPaths();
    
	// Draw final contour in output image
	imgContour = imgOutput;
	pPGComb->DrawInImageRGB(imgContour, iZoom);

	// Write final output image to BMP file
	if (imgContour.Save("./output.bmp", CImage2D::FORMAT_BMP)==false)
	{
		cerr<<"ERROR: unable to write to image file output_final.bmp"<<endl;
		return -1;
	}

    const CArray2D<float> &arrayPotential = pPGComb->GetPotential();
    CImage2D imgPotential;

    imgPotential.CreateFromArray2DFloat(arrayPotential);
    imgPotential.Save("./potential.jpg", CImage2D::FORMAT_JPG);

    unsigned idx;
    for (idx=0; idx<pPGComb->GetNbVertices(); idx++)
    {
        const list<CRealPath2D> &listPaths = pPGComb->GetAdmissibleSet(idx);
        char strFilename[100];
        imgContour = imgOutput;

        list<CRealPath2D>::const_iterator itPath;
        for (itPath = listPaths.begin(); itPath!=listPaths.end(); itPath++)
            itPath->DrawInImage(imgContour, iZoom);

        if (pPGComb->iNbAdmissiblePathsMax>1)
            pPGComb->GetMedialCurve(idx).DrawInImagePixelCorners(imgContour, iZoom);

        sprintf(strFilename, "./adm_paths%02d.bmp", idx);
        imgContour.Save(strFilename, CImage2D::FORMAT_BMP);

        imgContour = imgOutput;
        pPGComb->GetChosenPathInAdmissibleSet(idx).DrawInImage(imgContour, iZoom);

        sprintf(strFilename, "./chosen_path%02d.bmp", idx);
        imgContour.Save(strFilename, CImage2D::FORMAT_BMP);
    }

    // Make mask of inner region
    CImage2D imgMask;
    pPGComb->MakeBinaryMask(imgMask);
    imgMask.Save("./mask.bmp", CImage2D::FORMAT_BMP);

	// Destroy contour
	delete pPGComb;

	return 0;
}
