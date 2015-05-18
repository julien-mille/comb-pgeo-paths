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
histogram.cpp

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

#include <algorithm>
#include "histogram.h"

CHistogram1D::CHistogram1D()
{
	fGaussianStdDeviation = 0.0f;
    iNbBins = 0;
    iPadding = 0;
	fSumWeights = 0;
}

void CHistogram1D::SetNbBins(int nbBins)
{
    iNbBins = max(1, nbBins);
    Init(iNbBins + 2*iPadding);
    Fill(0.0f);
	fSumWeights = 0;
}

void CHistogram1D::SetGaussianStdDeviation(float stdDeviation)
{
    fGaussianStdDeviation = stdDeviation;
    arrayGaussian.SetGaussianKernel(stdDeviation);

    iPadding = (int)(3.0f*fGaussianStdDeviation);
	Init(iNbBins + 2*iPadding);
    Fill(0.0f);
	fSumWeights = 0;
}

void CHistogram1D::SetRange(float fMin, float fMax)
{
    fRangeMin = fMin;
    fRangeMax = fMax;
}

void CHistogram1D::AddElement(float fIntensity, float fWeight)
{
	int iIntensityScaled;

    #ifdef ARRAYND_RUNTIME_CHECK
    if (fIntensity<fRangeMin || fIntensity>fRangeMax)
    {
        cerr<<"ERROR in CHistogram1D::AddElement(...): value "<<fIntensity<<" out of range ["<<fRangeMin<<","<<fRangeMax<<"]"<<endl;
        return;
    }
    #endif
    iIntensityScaled = iPadding + (int)((fIntensity - fRangeMin)/(fRangeMax - fRangeMin) * (float)iNbBins);
    pElements[iIntensityScaled] += fWeight;
    fSumWeights += fWeight;
}

void CHistogram1D::SubtractElement(float fIntensity, float fWeight)
{
	int iIntensityScaled;

    #ifdef ARRAYND_RUNTIME_CHECK
    if (fIntensity<fRangeMin || fIntensity>fRangeMax)
    {
        cerr<<"ERROR in CHistogram1D::SubtractElement(...): value "<<fIntensity<<" out of range ["<<fRangeMin<<","<<fRangeMax<<"]"<<endl;
        return;
    }
    #endif
    iIntensityScaled = iPadding + (int)((fIntensity - fRangeMin)/(fRangeMax - fRangeMin) * (float)iNbBins);
    pElements[iIntensityScaled] -= fWeight;
    fSumWeights -= fWeight;
}

void CHistogram1D::AddElementScaled(int iIntensityScaled, float fWeight)
{
    #ifdef ARRAYND_RUNTIME_CHECK
    if (iIntensityScaled<iPadding || iIntensityScaled>=this->iSize-iPadding)
    {
        cerr<<"ERROR in CHistogram1D::AddElementScaled(...): scaled value "<<iIntensityScaled<<" out of range ["<<iPadding<<","<<this->iSize-1-iPadding<<"]"<<endl;
        return;
    }
    #endif
	pElements[iIntensityScaled] += fWeight;
	fSumWeights += fWeight;
}

void CHistogram1D::SubtractElementScaled(int iIntensityScaled, float fWeight)
{
    #ifdef ARRAYND_RUNTIME_CHECK
    if (iIntensityScaled<iPadding || iIntensityScaled>=this->iSize-iPadding)
    {
        cerr<<"ERROR in CHistogram1D::SubtractElementScaled(...): scaled value "<<iIntensityScaled<<" out of range ["<<iPadding<<","<<this->iSize-1-iPadding<<"]"<<endl;
        return;
    }
    #endif
	pElements[iIntensityScaled] -= fWeight;
	fSumWeights -= fWeight;
}

CHistogram1D &CHistogram1D::operator +=(const CHistogram1D &histo)
{
    int i;
    float *pValue;
    const float *pValue2;

    if (histo.iNbBins!=iNbBins || histo.iPadding!=iPadding)
	{
        cerr<<"ERROR in CHistogram1D::operator +=(...): histograms have different sizes"<<endl;
        return *this;
	}

    if (histo.fRangeMin!=fRangeMin || histo.fRangeMax!=fRangeMax)
	{
        cerr<<"ERROR in CHistogram1D::operator +=(...): histograms have different ranges"<<endl;
        return *this;
	}

    pValue = GetBuffer() + iPadding;
    pValue2 = histo.GetBuffer() + iPadding;

    for (i=0; i<iNbBins; i++)
    {
        *pValue += *pValue2;
        pValue++;
        pValue2++;
    }

    fSumWeights += histo.fSumWeights;

    return *this;
}

CHistogram1D &CHistogram1D::operator -=(const CHistogram1D &histo)
{
    int i;
    float *pValue;
    const float *pValue2;

    if (histo.iNbBins!=iNbBins || histo.iPadding!=iPadding)
	{
        cerr<<"ERROR in CHistogram1D::operator -=(...): histograms have different sizes"<<endl;
        return *this;
	}

    if (histo.fRangeMin!=fRangeMin || histo.fRangeMax!=fRangeMax)
	{
        cerr<<"ERROR in CHistogram1D::operator -=(...): histograms have different ranges"<<endl;
        return *this;
	}

    pValue = GetBuffer() + iPadding;
    pValue2 = histo.GetBuffer() + iPadding;

    for (i=0; i<iNbBins; i++)
    {
        *pValue -= *pValue2;
        pValue++;
        pValue2++;
    }

    fSumWeights -= histo.fSumWeights;

    return *this;
}

CHistogram1D CHistogram1D::operator +(const CHistogram1D &histo) const
{
    if (histo.iNbBins!=iNbBins || histo.iPadding!=iPadding)
	{
        cerr<<"ERROR in CHistogram1D::operator +(...): histograms have different sizes"<<endl;
        return *this;
	}

    if (histo.fRangeMin!=fRangeMin || histo.fRangeMax!=fRangeMax)
	{
        cerr<<"ERROR in CHistogram1D::operator +(...): histograms have different ranges"<<endl;
        return *this;
	}

    CHistogram1D histoRes;
    int i;
    const float *pValue, *pValue2;
    float *pValueRes;

    histoRes.iNbBins = iNbBins;
    histoRes.iPadding = iPadding;
    histoRes.fGaussianStdDeviation = fGaussianStdDeviation;
    histoRes.arrayGaussian = arrayGaussian;
    histoRes.fRangeMin = fRangeMin;
    histoRes.fRangeMax = fRangeMax;
    histoRes.fSumWeights = fSumWeights + histo.fSumWeights;

    histoRes.Init(iNbBins + 2*iPadding);

    pValue = GetBuffer();
    pValue2 = histo.GetBuffer();
    pValueRes = histoRes.GetBuffer();

    for (i=0; i<iSize; i++)
    {
        *pValueRes = *pValue + *pValue2;
        pValue++;
        pValue2++;
        pValueRes++;
    }

    return histoRes;
}

CHistogram1D CHistogram1D::operator -(const CHistogram1D &histo) const
{
    if (histo.iNbBins!=iNbBins || histo.iPadding!=iPadding)
	{
        cerr<<"ERROR in CHistogram1D::operator -(...): histograms have different sizes"<<endl;
        return *this;
	}

    if (histo.fRangeMin!=fRangeMin || histo.fRangeMax!=fRangeMax)
	{
        cerr<<"ERROR in CHistogram1D::operator -(...): histograms have different ranges"<<endl;
        return *this;
	}

    CHistogram1D histoRes;
    int i;
    const float *pValue, *pValue2;
    float *pValueRes;

    histoRes.iNbBins = iNbBins;
    histoRes.iPadding = iPadding;
    histoRes.fGaussianStdDeviation = fGaussianStdDeviation;
    histoRes.arrayGaussian = arrayGaussian;
    histoRes.fRangeMin = fRangeMin;
    histoRes.fRangeMax = fRangeMax;
    histoRes.fSumWeights = fSumWeights - histo.fSumWeights;

    histoRes.Init(iNbBins + 2*iPadding);

    pValue = GetBuffer();
    pValue2 = histo.GetBuffer();
    pValueRes = histoRes.GetBuffer();

    for (i=0; i<iSize; i++)
    {
        *pValueRes = *pValue - *pValue2;
        pValue++;
        pValue2++;
        pValueRes++;
    }

    return histoRes;
}

float CHistogram1D::GetValue(float fIntensity) const
{
    float fIntensityScaled;

    #ifdef ARRAYND_RUNTIME_CHECK
    if (fIntensity<fRangeMin || fIntensity>fRangeMax)
    {
        cerr<<"ERROR in CHistogram1D::SubtractElement(...): value "<<fIntensity<<" out of range ["<<fRangeMin<<","<<fRangeMax<<"]"<<endl;
        return 0.0f;
    }
    #endif
    fIntensityScaled = (float)iPadding + (fIntensity - fRangeMin)/(fRangeMax - fRangeMin) * (float)iNbBins;

    return GetElementInterpolate(fIntensityScaled);
}

void CHistogram1D::MakeImageScaled(const CArray2D<float> &arrayImage, CArray2D<int> &arrayImageScaled) const
{
    const float *pImage;
    int *pImageScaled;
    int i, iSize;
    bool bValueOutOfRange;

    arrayImageScaled.Init(arrayImage.GetSize());

    iSize = arrayImage.GetWidth()*arrayImage.GetHeight();
    pImage = arrayImage.GetBuffer();
    pImageScaled = arrayImageScaled.GetBuffer();
    bValueOutOfRange = false;
    for (i=0; i<iSize; i++)
    {
        if (*pImage>=fRangeMin && *pImage<=fRangeMax)
            *pImageScaled = iPadding + (int)((*pImage - fRangeMin)/(fRangeMax - fRangeMin) * (float)iNbBins);
        else {
            if (bValueOutOfRange==false)
            {
                cout<<"WARNING in CHistogram1D::MakeImageScaled(...): there were values out of range. Truncated them."<<endl;
                bValueOutOfRange = true;
            }
            if (*pImage<fRangeMin)
                *pImageScaled = iPadding;
            else
                *pImageScaled = this->iSize-1-iPadding;
        }
        pImage++;
        pImageScaled++;
    }
}

void CHistogram1D::GaussianSmooth()
{
    CArray1D<float>::ConvolveNoBorder(arrayGaussian);
}

float CHistogram1D::BhattacharyyaCoefficient(const CHistogram1D &histo) const
{
	const float *pValue, *pValue2;
	float fProba, fProba2;
	float fBhatCoef;
	int i;

	if (histo.iNbBins!=iNbBins || histo.iPadding!=iPadding)
	{
        cerr<<"ERROR in CHistogram1D::BhattacharyyaCoefficient(...): histograms have different sizes"<<endl;
        return 0.0f;
	}

	fBhatCoef = 0.0f;

    pValue = pElements + iPadding;
    pValue2 = histo.pElements + iPadding;

    for (i=0; i<iNbBins; i++)
    {
        if (*pValue<=0.0f)
            fProba = 0.0f;
        else
            fProba = *pValue/fSumWeights;

        if (*pValue2<=0.0f)
            fProba2 = 0.0f;
        else
            fProba2 = *pValue2/histo.fSumWeights;

        fBhatCoef += sqrt(fProba*fProba2);

        pValue++;
        pValue2++;
    }

	return fBhatCoef;
}



CHistogram3D::CHistogram3D()
{
	fGaussianStdDeviation = 0.0f;
    iNbBins.Set(0, 0, 0);
    iPadding = 0;
	fSumWeights = 0;
}

void CHistogram3D::SetNbBins(const CTriplet<int> &nbBins)
{
    iNbBins = tripletMax(CTriplet<int>(1,1,1), nbBins);
    Init(iNbBins + CTriplet<int>(iPadding, iPadding, iPadding)*2);
    Fill(0.0f);
	fSumWeights = 0;
}

void CHistogram3D::SetGaussianStdDeviation(float stdDeviation)
{
    fGaussianStdDeviation = stdDeviation;
    arrayGaussian.SetGaussianKernel(stdDeviation);

    iPadding = (int)(3.0f*fGaussianStdDeviation);
	Init(iNbBins + CTriplet<int>(iPadding, iPadding, iPadding)*2);
    Fill(0.0f);
	fSumWeights = 0;
}

void CHistogram3D::SetRange(const CTriplet<float> &fMin, const CTriplet<float> &fMax)
{
    fRangeMin = fMin;
    fRangeMax = fMax;
}

void CHistogram3D::AddElement(const CTriplet<float> &fColor, float fWeight)
{
	CTriplet<int> iColorScaled;

    #ifdef ARRAYND_RUNTIME_CHECK
    if (!fColor.IsInRange(fRangeMin, fRangeMax))
    {
        cerr<<"ERROR in CHistogram3D::AddElement(...): value "<<fColor;
        cerr<<" out of range ["<<fRangeMin.x<<".."<<fRangeMax.x<<"]x["<<fRangeMin.y<<".."<<fRangeMax.y<<"]x["<<fRangeMin.z<<".."<<fRangeMax.z<<"]"<<endl;
        return;
    }
    #endif
    iColorScaled.x = iPadding + (int)((fColor.x - fRangeMin.x)/(fRangeMax.x - fRangeMin.x) * (float)iNbBins.x);
    iColorScaled.y = iPadding + (int)((fColor.y - fRangeMin.y)/(fRangeMax.y - fRangeMin.y) * (float)iNbBins.y);
    iColorScaled.z = iPadding + (int)((fColor.z - fRangeMin.z)/(fRangeMax.z - fRangeMin.z) * (float)iNbBins.z);

    Element(iColorScaled) += fWeight;
    fSumWeights += fWeight;
}

void CHistogram3D::SubtractElement(const CTriplet<float> &fColor, float fWeight)
{
	CTriplet<int> iColorScaled;

    #ifdef ARRAYND_RUNTIME_CHECK
    if (!fColor.IsInRange(fRangeMin, fRangeMax))
    {
        cerr<<"ERROR in CHistogram3D::SubtractElement(...): value "<<fColor;
        cerr<<" out of range ["<<fRangeMin.x<<".."<<fRangeMax.x<<"]x["<<fRangeMin.y<<".."<<fRangeMax.y<<"]x["<<fRangeMin.z<<".."<<fRangeMax.z<<"]"<<endl;
        return;
    }
    #endif
    iColorScaled.x = iPadding + (int)((fColor.x - fRangeMin.x)/(fRangeMax.x - fRangeMin.x) * (float)iNbBins.x);
    iColorScaled.y = iPadding + (int)((fColor.y - fRangeMin.y)/(fRangeMax.y - fRangeMin.y) * (float)iNbBins.y);
    iColorScaled.z = iPadding + (int)((fColor.z - fRangeMin.z)/(fRangeMax.z - fRangeMin.z) * (float)iNbBins.z);

    Element(iColorScaled) -= fWeight;
    fSumWeights -= fWeight;
}

// Add and subtract scaled integer elements
void CHistogram3D::AddElementScaled(const CTriplet<int> &iColorScaled, float fWeight)
{
    #ifdef ARRAYND_RUNTIME_CHECK
    CTriplet<int> iRangeMinScaled(iPadding, iPadding, iPadding);
    CTriplet<int> iRangeMaxScaled = GetSize()-CTriplet<int>(1,1,1)-iRangeMinScaled;
    if (!iColorScaled.IsInRange(iRangeMinScaled, iRangeMaxScaled))
    {
        cerr<<"ERROR in CHistogram3D::AddElementScaled(...): scaled value "<<iColorScaled;
        cerr<<" out of range ["<<iRangeMinScaled.x<<".."<<iRangeMaxScaled.x<<"]x["<<iRangeMinScaled.y<<".."<<iRangeMaxScaled.y<<"]x["<<iRangeMinScaled.z<<".."<<iRangeMaxScaled.z<<"]"<<endl;
        return;
    }
    #endif
	Element(iColorScaled) += fWeight;
	fSumWeights += fWeight;
}

void CHistogram3D::SubtractElementScaled(const CTriplet<int> &iColorScaled, float fWeight)
{
    #ifdef ARRAYND_RUNTIME_CHECK
    CTriplet<int> iRangeMinScaled(iPadding, iPadding, iPadding);
    CTriplet<int> iRangeMaxScaled = GetSize()-CTriplet<int>(1,1,1)-iRangeMinScaled;
    if (!iColorScaled.IsInRange(iRangeMinScaled, iRangeMaxScaled))
    {
        cerr<<"ERROR in CHistogram3D::SubtractElementScaled(...): scaled value "<<iColorScaled;
        cerr<<" out of range ["<<iRangeMinScaled.x<<".."<<iRangeMaxScaled.x<<"]x["<<iRangeMinScaled.y<<".."<<iRangeMaxScaled.y<<"]x["<<iRangeMinScaled.z<<".."<<iRangeMaxScaled.z<<"]"<<endl;
        return;
    }
    #endif
	Element(iColorScaled) -= fWeight;
	fSumWeights -= fWeight;
}

CHistogram3D &CHistogram3D::operator +=(const CHistogram3D &histo)
{
    CTriplet<int> piCurrent;
    float *pValue;
    const float *pValue2;

    if (histo.iNbBins!=iNbBins || histo.iPadding!=iPadding)
	{
        cerr<<"ERROR in CHistogram3D::operator +=(...): histograms have different sizes"<<endl;
        return *this;
	}

    if (histo.fRangeMin!=fRangeMin || histo.fRangeMax!=fRangeMax)
	{
        cerr<<"ERROR in CHistogram3D::operator +=(...): histograms have different ranges"<<endl;
        return *this;
	}

	for (piCurrent.z=0; piCurrent.z<iNbBins.z; piCurrent.z++)
    {
        for (piCurrent.y=0; piCurrent.y<iNbBins.y; piCurrent.y++)
        {
            pValue = GetBuffer() + GetOffset(iPadding, iPadding+piCurrent.y, iPadding+piCurrent.z);
            pValue2 = histo.GetBuffer() + histo.GetOffset(iPadding, iPadding+piCurrent.y, iPadding+piCurrent.z);

            for (piCurrent.x=0; piCurrent.x<iNbBins.x; piCurrent.x++)
            {
                *pValue += *pValue2;
                pValue++;
                pValue2++;
            }
        }
    }
    fSumWeights += histo.fSumWeights;

    return *this;
}

CHistogram3D &CHistogram3D::operator -=(const CHistogram3D &histo)
{
    CTriplet<int> piCurrent;
    float *pValue;
    const float *pValue2;

    if (histo.iNbBins!=iNbBins || histo.iPadding!=iPadding)
	{
        cerr<<"ERROR in CHistogram3D::operator -=(...): histograms have different sizes"<<endl;
        return *this;
	}

    if (histo.fRangeMin!=fRangeMin || histo.fRangeMax!=fRangeMax)
	{
        cerr<<"ERROR in CHistogram3D::operator -=(...): histograms have different ranges"<<endl;
        return *this;
	}

	for (piCurrent.z=0; piCurrent.z<iNbBins.z; piCurrent.z++)
    {
        for (piCurrent.y=0; piCurrent.y<iNbBins.y; piCurrent.y++)
        {
            pValue = GetBuffer() + GetOffset(iPadding, iPadding+piCurrent.y, iPadding+piCurrent.z);
            pValue2 = histo.GetBuffer() + histo.GetOffset(iPadding, iPadding+piCurrent.y, iPadding+piCurrent.z);

            for (piCurrent.x=0; piCurrent.x<iNbBins.x; piCurrent.x++)
            {
                *pValue -= *pValue2;
                pValue++;
                pValue2++;
            }
        }
    }
    fSumWeights -= histo.fSumWeights;

    return *this;
}

CHistogram3D CHistogram3D::operator +(const CHistogram3D &histo) const
{
    if (histo.iNbBins!=iNbBins || histo.iPadding!=iPadding)
	{
        cerr<<"ERROR in CHistogram3D::operator +(...): histograms have different sizes"<<endl;
        return *this;
	}

    if (histo.fRangeMin!=fRangeMin || histo.fRangeMax!=fRangeMax)
	{
        cerr<<"ERROR in CHistogram3D::operator +(...): histograms have different ranges"<<endl;
        return *this;
	}

    CHistogram3D histoRes;
    int i, iSize;
    const float *pValue, *pValue2;
    float *pValueRes;

    histoRes.iNbBins = iNbBins;
    histoRes.iPadding = iPadding;
    histoRes.fGaussianStdDeviation = fGaussianStdDeviation;
    histoRes.arrayGaussian = arrayGaussian;
    histoRes.fRangeMin = fRangeMin;
    histoRes.fRangeMax = fRangeMax;
    histoRes.fSumWeights = fSumWeights + histo.fSumWeights;

    histoRes.Init(iNbBins + CTriplet<int>(iPadding, iPadding, iPadding)*2);
    iSize = CArray1D<float>::GetSize();

    pValue = GetBuffer();
    pValue2 = histo.GetBuffer();
    pValueRes = histoRes.GetBuffer();

    for (i=0; i<iSize; i++)
    {
        *pValueRes = *pValue + *pValue2;
        pValue++;
        pValue2++;
        pValueRes++;
    }

    return histoRes;
}

CHistogram3D CHistogram3D::operator -(const CHistogram3D &histo) const
{
    if (histo.iNbBins!=iNbBins || histo.iPadding!=iPadding)
	{
        cerr<<"ERROR in CHistogram3D::operator -(...): histograms have different sizes"<<endl;
        return *this;
	}

    if (histo.fRangeMin!=fRangeMin || histo.fRangeMax!=fRangeMax)
	{
        cerr<<"ERROR in CHistogram3D::operator -(...): histograms have different ranges"<<endl;
        return *this;
	}

    CHistogram3D histoRes;
    int i, iSize;
    const float *pValue, *pValue2;
    float *pValueRes;

    histoRes.iNbBins = iNbBins;
    histoRes.iPadding = iPadding;
    histoRes.fGaussianStdDeviation = fGaussianStdDeviation;
    histoRes.arrayGaussian = arrayGaussian;
    histoRes.fRangeMin = fRangeMin;
    histoRes.fRangeMax = fRangeMax;
    histoRes.fSumWeights = fSumWeights - histo.fSumWeights;

    histoRes.Init(iNbBins + CTriplet<int>(iPadding, iPadding, iPadding)*2);
    iSize = CArray1D<float>::GetSize();

    pValue = GetBuffer();
    pValue2 = histo.GetBuffer();
    pValueRes = histoRes.GetBuffer();

    for (i=0; i<iSize; i++)
    {
        *pValueRes = *pValue - *pValue2;
        pValue++;
        pValue2++;
        pValueRes++;
    }

    return histoRes;
}

float CHistogram3D::GetValue(const CTriplet<float> &fColor) const
{
	CTriplet<float> fColorScaled;

    #ifdef ARRAYND_RUNTIME_CHECK
    if (!fColor.IsInRange(fRangeMin, fRangeMax))
    {
        cerr<<"ERROR in CHistogram3D::GetValue(...): value "<<fColor;
        cerr<<" out of range ["<<fRangeMin.x<<".."<<fRangeMax.x<<"]x["<<fRangeMin.y<<".."<<fRangeMax.y<<"]x["<<fRangeMin.z<<".."<<fRangeMax.z<<"]"<<endl;
        return 0.0f;
    }
    #endif
    fColorScaled.x = (float)iPadding + (fColor.x - fRangeMin.x)/(fRangeMax.x - fRangeMin.x) * (float)iNbBins.x;
    fColorScaled.y = (float)iPadding + (fColor.y - fRangeMin.y)/(fRangeMax.y - fRangeMin.y) * (float)iNbBins.y;
    fColorScaled.z = (float)iPadding + (fColor.z - fRangeMin.z)/(fRangeMax.z - fRangeMin.z) * (float)iNbBins.z;

    return GetElementInterpolate(fColorScaled);
}

void CHistogram3D::MakeImageScaled(const CArray2D<CTriplet<float> > &arrayImage,
    CArray2D<CTriplet<int> > &arrayImageScaled) const
{
    const CTriplet<float> *pImage;
    CTriplet<int> *pImageScaled;
    int i, iImageSize;
    bool bValueOutOfRange;

    CTriplet<int> iRangeMinScaled(iPadding, iPadding, iPadding);
    CTriplet<int> iRangeMaxScaled = GetSize()-CTriplet<int>(1,1,1)-iRangeMinScaled;

    arrayImageScaled.Init(arrayImage.GetSize());

    iImageSize = arrayImage.GetWidth()*arrayImage.GetHeight();
    pImage = arrayImage.GetBuffer();
    pImageScaled = arrayImageScaled.GetBuffer();
    bValueOutOfRange = false;
    for (i=0; i<iImageSize; i++)
    {
        if (pImage->x==fRangeMax.x)
            pImageScaled->x = iPadding + iNbBins.x-1;
        else
            pImageScaled->x = iPadding + (int)((pImage->x - fRangeMin.x)/(fRangeMax.x - fRangeMin.x) * (float)iNbBins.x);

        if (pImage->y==fRangeMax.y)
            pImageScaled->y = iPadding + iNbBins.y-1;
        else
            pImageScaled->y = iPadding + (int)((pImage->y - fRangeMin.y)/(fRangeMax.y - fRangeMin.y) * (float)iNbBins.y);

        if (pImage->z==fRangeMax.z)
            pImageScaled->z = iPadding + iNbBins.z-1;
        else
            pImageScaled->z = iPadding + (int)((pImage->z - fRangeMin.z)/(fRangeMax.z - fRangeMin.z) * (float)iNbBins.z);

        if (!pImageScaled->IsInRange(iRangeMinScaled, iRangeMaxScaled))
        {
            if (bValueOutOfRange==false)
            {
                cout<<"WARNING in CHistogram3D::MakeImageScaled(...): there were values out of range. Truncated them."<<endl;
                bValueOutOfRange = true;
            }
            pImageScaled->LimitWithinRange(iRangeMinScaled, iRangeMaxScaled);
        }

        pImage++;
        pImageScaled++;
    }
}

void CHistogram3D::GaussianSmooth()
{
    ConvolveInPlaceX(arrayGaussian);
    ConvolveInPlaceY(arrayGaussian);
    ConvolveInPlaceZ(arrayGaussian);
}

float CHistogram3D::BhattacharyyaCoefficient(const CHistogram3D &histo) const
{
	const float *pValue, *pValue2;
	float fProba, fProba2;
	CTriplet<int> piCurrent;
	float fBhatCoef;

	if (histo.iNbBins!=iNbBins || histo.iPadding!=iPadding)
	{
        cerr<<"ERROR in CHistogram3D::BhattacharyyaCoefficient(...): histograms have different sizes"<<endl;
        return 0.0f;
	}

    if (histo.fRangeMin!=fRangeMin || histo.fRangeMax!=fRangeMax)
	{
        cerr<<"ERROR in CHistogram3D::BhattacharyyaCoefficient(...): histograms have different ranges"<<endl;
        return 0.0f;
	}

	fBhatCoef = 0.0f;

    for (piCurrent.z=0; piCurrent.z<iNbBins.z; piCurrent.z++)
    {
        for (piCurrent.y=0; piCurrent.y<iNbBins.y; piCurrent.y++)
        {
            pValue = GetBuffer() + GetOffset(iPadding, iPadding+piCurrent.y, iPadding+piCurrent.z);
            pValue2 = histo.GetBuffer() + histo.GetOffset(iPadding, iPadding+piCurrent.y, iPadding+piCurrent.z);

            for (piCurrent.x=0; piCurrent.x<iNbBins.x; piCurrent.x++)
            {
                if (*pValue<=0.0f)
                    fProba = 0.0f;
                else
                    fProba = *pValue/fSumWeights;

                if (*pValue2<=0.0f)
                    fProba2 = 0.0f;
                else
                    fProba2 = *pValue2/histo.fSumWeights;

                fBhatCoef += sqrt(fProba*fProba2);

                pValue++;
                pValue2++;
            }
        }
    }

	return fBhatCoef;
}
