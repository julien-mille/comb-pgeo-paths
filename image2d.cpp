#include "image2d.h"
#include <iostream>
#include <vector>
#include <stdio.h>
#include <string.h>

#ifdef IMAGE2D_SUPPORT_PNG
extern "C" {
#include <png.h>
}
#endif

#ifdef IMAGE2D_SUPPORT_JPG
extern "C" {
#include <jpeglib.h>
}
#include <setjmp.h>
#endif

#ifdef _MSC_VER
#define strcasecmp _stricmp
#endif

using namespace std;

CArray1D<unsigned char> CImage2D::vectBitsBitmapInfo;

#if !defined(_WINDOWS_) && !defined(_WINDOWS_H)
typedef struct {
	unsigned char rgbBlue, rgbGreen, rgbRed, rgbReserved;
} RGBQUAD;

typedef unsigned int IMAGE2D_DWORD;
typedef int IMAGE2D_LONG;
typedef unsigned short IMAGE2D_WORD;

typedef struct {
	IMAGE2D_DWORD   bcSize;
	IMAGE2D_WORD    bcWidth;
	IMAGE2D_WORD    bcHeight;
	IMAGE2D_WORD    bcPlanes;
	IMAGE2D_WORD    bcBitCount;
} BITMAPCOREHEADER;

typedef struct {
	IMAGE2D_DWORD      biSize;
	IMAGE2D_LONG       biWidth;
	IMAGE2D_LONG       biHeight;
	IMAGE2D_WORD       biPlanes;
	IMAGE2D_WORD       biBitCount;
	IMAGE2D_DWORD      biCompression;
	IMAGE2D_DWORD      biSizeImage;
	IMAGE2D_LONG       biXPelsPerMeter;
	IMAGE2D_LONG       biYPelsPerMeter;
	IMAGE2D_DWORD      biClrUsed;
	IMAGE2D_DWORD      biClrImportant;
} BITMAPINFOHEADER;

// Save default alignment rule and set it to 2 bytes
// so that sizeof(BITMAPFILEHADER) will return actual size of structure (14 bytes)
#pragma pack(push, r1, 2)
typedef struct {
        IMAGE2D_WORD    bfType;
        IMAGE2D_DWORD   bfSize;
        IMAGE2D_WORD    bfReserved1;
        IMAGE2D_WORD    bfReserved2;
        IMAGE2D_DWORD   bfOffBits;
} BITMAPFILEHEADER;
#pragma pack(pop, r1) // Restore default alignment rule

// Values for BITMAPINFOHEADER's biCompression field
#define BI_RGB        0L
#define BI_RLE8       1L
#define BI_RLE4       2L
#define BI_BITFIELDS  3L
#define BI_JPEG       4L
#define BI_PNG        5L
#endif // !defined(_WINDOWS_) && !defined(_WINDOWS_H)

#ifdef IMAGE2D_RUNTIME_CHECK
// Pixel returned in case of invalid access with function Pixel()
CImage2DPixel CImage2D::pixelError;
#endif

void CImage2D::StartImage2D()
{
	RGBQUAD *pPalette;
	unsigned char iColor;

	CImage2D::vectBitsBitmapInfo.Init(sizeof(BITMAPINFOHEADER)+sizeof(RGBQUAD)*256);

	pPalette = (RGBQUAD *)(CImage2D::vectBitsBitmapInfo.GetBuffer()+sizeof(BITMAPINFOHEADER));
	for (iColor=0;iColor<255;iColor++)
	{
		pPalette->rgbRed = iColor;
		pPalette->rgbGreen = iColor;
		pPalette->rgbBlue = iColor;
		pPalette->rgbReserved = 0;
		pPalette++;
	}
	pPalette->rgbRed = iColor;
	pPalette->rgbGreen = iColor;
	pPalette->rgbBlue = iColor;
	pPalette->rgbReserved = 0;
}

CImage2D::CImage2D(const CImage2D &image):CArray1D<unsigned char>(image)
{
	iSize = image.iSize;
	iWidth = image.iWidth;
	iHeight = image.iHeight;
	iBytesPerLine = image.iBytesPerLine;
	iBitsPerPixel = image.iBitsPerPixel;

	pElements = new CImage2DPixel[iSize];
	memcpy(pElements, image.pElements, iSize*sizeof(CImage2DPixel));
}


bool CImage2D::Create(int width, int height, int bitsPerPixel)
{
	int iBytesPerLineTemp;

	// Test if number of bits per pixel is valid
	if (bitsPerPixel!=8 && bitsPerPixel!=24 && bitsPerPixel!=32)
	{
		cerr<<"ERROR in CImage2D::Create(int,int,int): number of bits/pixel is not valid"<<endl;
		Empty();
		return false;
	}

	iWidth = width;
	iHeight = height;
	iBitsPerPixel = bitsPerPixel;

	iBytesPerLineTemp = iWidth*iBitsPerPixel/8;

	if (iBytesPerLineTemp%4!=0)
		iBytesPerLine = iBytesPerLineTemp+(4-iBytesPerLineTemp%4);
	else iBytesPerLine = iBytesPerLineTemp;

	if (CArray1D<CImage2DPixel>::Init(iBytesPerLine*iHeight)==false)
	{
		iWidth = 0;
		iHeight = 0;
		iBitsPerPixel = 0;
		return false;
	}
	return true;
}

bool CImage2D::Load(const char *strFileName)
{
	char strExtension[20];
	strcpy(strExtension, strFileName + strlen(strFileName) - 4);

	if (strcasecmp(strExtension, ".bmp")==0)
	{
		FILE *pFile;

		pFile = fopen(strFileName, "rb");
		if (pFile==NULL)
		{
			cerr<<"ERROR in CImage2D::Load(...): cannot open file "<<strFileName<<endl;
			Empty();
			return false;
		}

		if (DecodeBMP(pFile)==false)
		{
			fclose(pFile);
			Empty();
			return false;
		}

		fclose(pFile);
		return true;
	}
#ifdef IMAGE2D_SUPPORT_JPG
	if (strcasecmp(strExtension, ".jpg")==0)
	{
		FILE *pFile;

		pFile = fopen(strFileName, "rb");
		if (pFile==NULL)
		{
			cerr<<"ERROR in CImage2D::Load(...): cannot open file "<<strFileName<<endl;
			Empty();
			return false;
		}

		if (DecodeJPG(pFile)==false)
		{
			fclose(pFile);
			Empty();
			return false;
		}

		fclose(pFile);
		return true;
	}
#endif
#ifdef IMAGE2D_SUPPORT_PNG
	if (strcasecmp(strExtension, ".png")==0)
	{
		FILE *pFile;

		pFile = fopen(strFileName, "rb");
		if (pFile==NULL)
		{
			cerr<<"ERROR in CImage2D::Load(...): cannot open file "<<strFileName<<endl;
			Empty();
			return false;
		}
		if (DecodePNG(pFile)==false)
		{
			fclose(pFile);
			Empty();
			return false;
		}

		fclose(pFile);
		return true;
	}
#endif

	// Image format is not supported
	cerr<<"ERROR in CImage2D::Load(...): image format unrecognized or not supported"<<endl;
	Empty();
	return false;
}

bool CImage2D::Save(const char *strFileName, typeFormat format) const
{
	if (format==FORMAT_BMP)
	{
		FILE *pFile;
		bool bReturn;

		pFile = fopen(strFileName, "wb");
		if (pFile==NULL)
		{
			cerr<<"ERROR in CImage2D::Save(...): cannot open file "<<strFileName<<" for writing"<<endl;
			return false;
		}

		bReturn = EncodeBMP(pFile);

		fclose(pFile);
		return bReturn;
	}

#ifdef IMAGE2D_SUPPORT_JPG
	if (format==FORMAT_JPG)
	{
		FILE *pFile;
		bool bReturn;

		pFile = fopen(strFileName, "wb");
		if (pFile==NULL)
		{
			cerr<<"ERROR in CImage2D::Save(...): cannot open file "<<strFileName<<" for writing"<<endl;
			return false;
		}
		bReturn = EncodeJPG(pFile);

		fclose(pFile);
		return bReturn;
	}
#endif
#ifdef IMAGE2D_SUPPORT_PNG
	if (format==FORMAT_PNG)
	{
		FILE *pFile;
		bool bReturn;

		pFile = fopen(strFileName, "wb");
		if (pFile==NULL)
		{
			cerr<<"ERROR in CImage2D::Save(...): cannot open file "<<strFileName<<" for writing"<<endl;
			return false;
		}
		bReturn = EncodePNG(pFile);

		fclose(pFile);
		return bReturn;
	}
#endif

	// Image format is not supported
	cerr<<"ERROR in CImage2D::Save(...): image format unrecognized or not supported"<<endl;
	return false;
}

void CImage2D::Empty()
{
	iWidth = iHeight = iBytesPerLine = iBitsPerPixel = 0;
	CArray1D<CImage2DPixel>::Empty();
}

CImage2D &CImage2D::operator =(const CImage2D &image)
{
	CArray1D<CImage2DPixel>::operator =(image);

	iWidth = image.iWidth;
	iHeight = image.iHeight;
	iBytesPerLine = image.iBytesPerLine;
	iBitsPerPixel = image.iBitsPerPixel;

	return *this;
}

CImage2D CImage2D::Resize(int newWidth, int newHeight) const
{
	CImage2D dest;
	int x, y, iOffsetX, iOffsetEndRowDest;
	CImage2DPixel *pRow, *pBitsDest;
	float fScaleX, fScaleY;

	#ifdef IMAGE2D_RUNTIME_CHECK
	if (iBitsPerPixel!=8)
	{
		cerr<<"ERROR in CImage2D::Resize(...): image is not 8-bit. Call ResizeRGB() instead ?"<<endl;
		return dest;
	}
	#endif

	if (dest.Create(newWidth, newHeight, 8)==false)
	{
		cerr<<"ERROR in CImage2D::Resize(...): cannot create image"<<endl;
		return dest;
	}

	iOffsetEndRowDest = dest.iBytesPerLine - dest.iWidth;

	fScaleX = ((float)iWidth)/((float)newWidth);
	fScaleY = ((float)iHeight)/((float)newHeight);

	pBitsDest = dest.pElements;

	for (y=0; y<newHeight; y++)
	{
		pRow = pElements + (int)(((float)y)*fScaleY)*iBytesPerLine;
		for (x=0; x<newWidth; x++)
		{
			iOffsetX = (int)(((float)x)*fScaleX);
			*pBitsDest = pRow[iOffsetX];
			pBitsDest++;
		}
		pBitsDest += iOffsetEndRowDest;
	}
	return dest;
}

CImage2D CImage2D::ResizeRGB(int newWidth, int newHeight) const
{
	CImage2D dest;
	int x, y, iOffsetX, iOffsetEndRowDest, iBytesPerPixel;
	CImage2DPixel *pRow, *pBitsDest;
	float fScaleX, fScaleY;

	#ifdef IMAGE2D_RUNTIME_CHECK
	if (iBitsPerPixel!=24 && iBitsPerPixel!=32)
	{
		cerr<<"ERROR in CImage2D::ResizeRGB(...): image is not RGB. Call Resize() instead ?"<<endl;
		return dest;
	}
	#endif

	if (dest.Create(newWidth, newHeight, iBitsPerPixel)==false)
	{
		cerr<<"ERROR in CImage2D::ResizeRGB(...): cannot create image"<<endl;
		return dest;
	}

	iBytesPerPixel = iBitsPerPixel/8;
	iOffsetEndRowDest = dest.iBytesPerLine - dest.iWidth*iBytesPerPixel;

	fScaleX = ((float)iWidth)/((float)newWidth);
	fScaleY = ((float)iHeight)/((float)newHeight);

	pBitsDest = dest.pElements;

	for (y=0; y<newHeight; y++)
	{
		pRow = pElements + (int)(((float)y)*fScaleY)*iBytesPerLine;
		for (x=0; x<newWidth; x++)
		{
			iOffsetX = (int)(((float)x)*fScaleX)*iBytesPerPixel;
			*((CImage2DByteRGBPixel *)pBitsDest) = *((CImage2DByteRGBPixel *)(pRow+iOffsetX));
			pBitsDest += iBytesPerPixel;
		}
		pBitsDest += iOffsetEndRowDest;
	}
	return dest;
}

CImage2D CImage2D::GrayScale() const
{
	CImage2D dest;

	if (iBitsPerPixel==8)
	{
		cout<<"WARNING in CImage2D::Grayscale(): image is already 8-bit"<<endl;
		return *this;
	}

	if (iBitsPerPixel!=24 && iBitsPerPixel!=32)
	{
		cerr<<"ERROR in CImage2D::Grayscale(): nb bits per pixel not supported"<<endl;
		return dest;
	}

	int x, y, iOffsetEndRow, iOffsetEndRowDest, iBytesPerPixel;
	CImage2DPixel *pBitsDest, *pBits;

	if (!dest.Create(iWidth, iHeight, 8))
	{
		cerr<<"ERROR in CImage2D::Grayscale(): cannot create image"<<endl;
		return dest;
	}

	pBits = pElements;
	pBitsDest = dest.pElements;

	iBytesPerPixel = iBitsPerPixel/8;
	iOffsetEndRow = iBytesPerLine - iWidth*iBytesPerPixel;
	iOffsetEndRowDest = dest.iBytesPerLine - dest.iWidth;

	for (y=0;y<iHeight;y++)
	{
		for (x=0;x<iWidth;x++)
		{
			// Intensity is taken as the average of RGB components
			*pBitsDest = (unsigned char)(((int)pBits[0]+(int)pBits[1]+(int)pBits[2])/3);

			pBits+=iBytesPerPixel;
			pBitsDest++;
		}
		pBits+=iOffsetEndRow;
		pBitsDest+=iOffsetEndRowDest;
	}

	return dest;
}

bool CImage2D::IsGrayScaleRGB() const
{
	bool bIsGrayScale = false;

	if (iBitsPerPixel==24 || iBitsPerPixel==32)
	{
		int x, y, iOffsetEndRow, iBytesPerPixel;
		CImage2DPixel *pBits;

		pBits = pElements;

		iBytesPerPixel = iBitsPerPixel/8;
		iOffsetEndRow = iBytesPerLine - iWidth*iBytesPerPixel;

		bIsGrayScale = true;
		for (y=0;y<iHeight && bIsGrayScale==true;y++)
		{
			for (x=0;x<iWidth && bIsGrayScale==true;x++)
			{
				if (pBits[0]!=pBits[1] || pBits[1]!=pBits[2])
					bIsGrayScale = false;
				pBits+=iBytesPerPixel;
			}
			pBits+=iOffsetEndRow;
		}
	}
	return bIsGrayScale;
}

void CImage2D::Convert8bitsto24bits(CImage2D &dest24bits) const
{
	int x, y, width, height, dec8bits, dec24bits;
	CImage2DPixel *pPixelsSrc, *pPixels24bitsDest;

	width = iWidth;
	height = iHeight;

	if (dest24bits.GetWidth()!=width || dest24bits.GetHeight()!=height || dest24bits.GetBitsPerPixel()!=24)
		dest24bits.Create(width, height, 24);

	dec8bits = this->GetBytesPerLine()-width;
	dec24bits = dest24bits.GetBytesPerLine()-3*width;

	pPixelsSrc = this->GetBits();
	pPixels24bitsDest = dest24bits.GetBits();

	for (y=0;y<height;y++)
	{
		for (x=0;x<width;x++)
		{
			pPixels24bitsDest[0] = pPixels24bitsDest[1] = pPixels24bitsDest[2] = *pPixelsSrc;
			pPixelsSrc++;
			pPixels24bitsDest+=3;
		}
		pPixelsSrc        += dec8bits;
		pPixels24bitsDest += dec24bits;
	}
}

void CImage2D::TransformWithLUT(const CArray1D<unsigned char> &vectLUT)
{
	CImage2DPixel *pBits = pElements;
	int x, y, iOffsetEndRow;

	iOffsetEndRow = iBytesPerLine-iWidth;

	for (y=0;y<iHeight;y++)
	{
		for (x=0;x<iWidth;x++)
		{
			*pBits = vectLUT[(int)(*pBits)];
			pBits++;
		}
		pBits+=iOffsetEndRow;
	}
}

void CImage2D::TransformWithFunctionRGB(typeByteRGBPixelTransformFunction transformFunction)
{
	CImage2DPixel *pBits;
	int x, y, iOffsetEndRow, iBytesPerPixel;

	pBits = pElements;
	iBytesPerPixel = iBitsPerPixel/8;
	iOffsetEndRow = iBytesPerLine - iWidth*iBytesPerPixel;

	for (y=0; y<iHeight; y++)
	{
		for (x=0; x<iWidth; x++)
		{
			*((CImage2DByteRGBPixel *)pBits) = (*transformFunction)(*((CImage2DByteRGBPixel *)pBits));
			pBits += iBytesPerPixel;
		}
		pBits += iOffsetEndRow;
	}
}

CImage2D CImage2D::AbsDiff(const CImage2D &imgArg) const
{
	CImage2D imgDest;
	CImage2DPixel *pBits, *pBitsArg, *pBitsDest;
	int x, y, iOffsetEndRow;

    imgDest.Create(iWidth, iHeight, 8);

	iOffsetEndRow = iBytesPerLine-iWidth;
    pBits = pElements;
    pBitsArg = imgArg.pElements;
    pBitsDest = imgDest.pElements;

	for (y=0;y<iHeight;y++)
	{
		for (x=0;x<iWidth;x++)
		{
			*pBitsDest = (unsigned char)abs((int)(*pBits)-(int)(*pBitsArg));
			pBits++;
			pBitsArg++;
			pBitsDest++;
		}
		pBits+=iOffsetEndRow;
		pBitsArg+=iOffsetEndRow;
		pBitsDest+=iOffsetEndRow;
	}

	return imgDest;
}

CImage2D CImage2D::Min(const CImage2D &imgArg) const
{
	CImage2D imgDest;
	CImage2DPixel *pBits, *pBitsArg, *pBitsDest;
	int x, y, iOffsetEndRow;

    imgDest.Create(iWidth, iHeight, 8);

	iOffsetEndRow = iBytesPerLine-iWidth;
    pBits = pElements;
    pBitsArg = imgArg.pElements;
    pBitsDest = imgDest.pElements;

	for (y=0;y<iHeight;y++)
	{
		for (x=0;x<iWidth;x++)
		{
			*pBitsDest = min(*pBits, *pBitsArg);
			pBits++;
			pBitsArg++;
			pBitsDest++;
		}
		pBits+=iOffsetEndRow;
		pBitsArg+=iOffsetEndRow;
		pBitsDest+=iOffsetEndRow;
	}

	return imgDest;
}

CImage2D CImage2D::Max(const CImage2D &imgArg) const
{
	CImage2D imgDest;
	CImage2DPixel *pBits, *pBitsArg, *pBitsDest;
	int x, y, iOffsetEndRow;

    imgDest.Create(iWidth, iHeight, 8);

	iOffsetEndRow = iBytesPerLine-iWidth;
    pBits = pElements;
    pBitsArg = imgArg.pElements;
    pBitsDest = imgDest.pElements;

	for (y=0;y<iHeight;y++)
	{
		for (x=0;x<iWidth;x++)
		{
			*pBitsDest = max(*pBits, *pBitsArg);
			pBits++;
			pBitsArg++;
			pBitsDest++;
		}
		pBits+=iOffsetEndRow;
		pBitsArg+=iOffsetEndRow;
		pBitsDest+=iOffsetEndRow;
	}

	return imgDest;
}


void CImage2D::FloodFill(const CCouple<int> &piSeed, CImage2DPixel col)
{
	CImage2DPixel col_seed;
	CArray2D<unsigned char> arrayMarked;
	unsigned char *pMarked, *pMarkedNeighbor;
	unsigned char *pPixel, *pPixelNeighbor;
	int x, y, v, offset;

	CCouple<int> piCurrent, piNeighbor, piMin, piMax;
	vector<CCouple<int> > listPoints;
	unsigned int iPointCurrent;

	listPoints.reserve(iWidth*iHeight);

	piMin.Set(0, 0);
	piMax.Set(iWidth-1, iHeight-1);

	// Voisinage 2D (4-connexe)
	CArray1D<CCouple<int> > voisinage2D(4);
	voisinage2D[0].Set(-1,0);
	voisinage2D[1].Set(1,0);
	voisinage2D[2].Set(0,-1);
	voisinage2D[3].Set(0,1);

	// Offset correspondant au voisinage 2D (4-connexe)
	CArray1D<int> offsets_voisinage2D(4);
	offsets_voisinage2D[0] = -1;
	offsets_voisinage2D[1] = 1;
	offsets_voisinage2D[2] = -iWidth;
	offsets_voisinage2D[3] = iWidth;

	CArray1D<int> offsets_voisinage_image2D(4);
	offsets_voisinage_image2D[0] = -1;
	offsets_voisinage_image2D[1] = 1;
	offsets_voisinage_image2D[2] = -iBytesPerLine;
	offsets_voisinage_image2D[3] = iBytesPerLine;

	arrayMarked.Init(iWidth, iHeight);
	arrayMarked.Fill(2);

	// Erase bit 1 on borders : not in safety area
	for (y=0; y<iHeight; y++)
	{
		arrayMarked.Element(0,y) &= 0xFD;
		arrayMarked.Element(iWidth-1,y) &= 0xFD;
	}
	for (x=0; x<iWidth; x++)
	{
		arrayMarked.Element(x,0) &= 0xFD;
		arrayMarked.Element(x,iHeight-1) &= 0xFD;
	}

	// Colorie le germe initial
	offset = iBytesPerLine*piSeed.y+piSeed.x;
	col_seed = pElements[offset];
	pElements[offset] = col;
	arrayMarked.Element(piSeed) |= 0x01;

	listPoints.push_back(piSeed);
	iPointCurrent = 0;
	while (iPointCurrent<listPoints.size())
	{
		piCurrent = listPoints[iPointCurrent++];
        pPixel = pElements + iBytesPerLine*piCurrent.y + piCurrent.x;
        pMarked = arrayMarked.GetBuffer() + iWidth*piCurrent.y + piCurrent.x;

        if (((*pMarked)&0x02)!=0)
        {
            // Safety area -> no need to check neighbors
            for (v=0;v<voisinage2D.GetSize();v++)
            {
                pMarkedNeighbor = pMarked + offsets_voisinage2D[v];
                pPixelNeighbor = pPixel + offsets_voisinage_image2D[v];
                if (((*pMarkedNeighbor)&0x01)==0 && *pPixelNeighbor==col_seed) // [offsetv]==0)
                {
                    *pMarkedNeighbor |= 0x01;
                    *pPixelNeighbor = col;
                    listPoints.push_back(piCurrent+voisinage2D[v]);
                }
            }
        }
        else {
            for (v=0;v<voisinage2D.GetSize();v++)
            {
                piNeighbor = piCurrent+voisinage2D[v];
                if (piNeighbor.IsInRange(piMin, piMax))
                {
                    pMarkedNeighbor = pMarked + offsets_voisinage2D[v];
                    pPixelNeighbor = pPixel + offsets_voisinage_image2D[v];
                    if (((*pMarkedNeighbor)&0x01)==0 && *pPixelNeighbor==col_seed)
                    {
                        *pMarkedNeighbor |= 0x01;
                        *pPixelNeighbor = col;
                        listPoints.push_back(piNeighbor);
                    }
                }
            }
        }
	}
}

bool CImage2D::CreateFromArray2DFloat(const CArray2D<float> &array2D)
{
	int x, y, iOffsetEndRow;
	CImage2DPixel *pBits;
	const float *pFloat;

	if (iWidth!=array2D.GetWidth() || iHeight!=array2D.GetHeight() || iBitsPerPixel!=8)
	{
		if (Create(array2D.GetWidth(), array2D.GetHeight(), 8)==false)
			return false;
	}

	iOffsetEndRow = iBytesPerLine - iWidth;
	pFloat = array2D.GetBuffer();
	pBits = pElements;

	for (y=0;y<iHeight;y++)
	{
		for (x=0;x<iWidth;x++)
		{
			if (*pFloat<0.0f)
				*pBits = 0;
			else if (*pFloat>1.0f)
				*pBits = 255;
			else
				*pBits = (CImage2DPixel)(*pFloat*255.0f);
			pBits++;
			pFloat++;
		}
		pBits += iOffsetEndRow;
	}
	return true;
}

bool CImage2D::CreateFromArray2DTripletFloatRGB(CArray2D<CTriplet<float> > &array2D)
{
	int x, y, iOffsetEndRow;
	CImage2DPixel *pBits;
	CTriplet<float> *pTripletFloat;

	if (iWidth!=array2D.GetWidth() || iHeight!=array2D.GetHeight() || iBitsPerPixel!=24)
	{
		if (Create(array2D.GetWidth(), array2D.GetHeight(), 24)==false)
			return false;
	}

	pTripletFloat = array2D.GetBuffer();
	iOffsetEndRow = iBytesPerLine - iWidth*3;
	pBits = pElements;

	for (y=0;y<iHeight;y++)
	{
		for (x=0;x<iWidth;x++)
		{
			if (pTripletFloat->x<0.0f) pBits[2] = 0;
			else if (pTripletFloat->x>1.0f) pBits[2] = 255;
			else pBits[2] = (CImage2DPixel)(pTripletFloat->x*255.0f);

			if (pTripletFloat->y<0.0f) pBits[1] = 0;
			else if (pTripletFloat->y>1.0f) pBits[1] = 255;
			else pBits[1] = (CImage2DPixel)(pTripletFloat->y*255.0f);

			if (pTripletFloat->z<0.0f) pBits[0] = 0;
			else if (pTripletFloat->z>1.0f) pBits[0] = 255;
			else pBits[0] = (CImage2DPixel)(pTripletFloat->z*255.0f);

			pBits+=3;
			pTripletFloat++;
		}
		pBits += iOffsetEndRow;
	}
	return true;
}

bool CImage2D::ConvertToArray2DFloat(CArray2D<float> &array2D) const
{
	int x, y, iOffsetEndRow;
	CImage2DPixel *pBits;
	float *pFloat;

	#ifdef IMAGE2D_RUNTIME_CHECK
	if (iBitsPerPixel!=8)
	{
		cerr<<"ERROR in CImage2D::ConvertToArray2DFloat(...): image is not 8-bit"<<endl;
		return false;
	}
	#endif

	if (iWidth!=array2D.GetWidth() || iHeight!=array2D.GetHeight())
	{
		if (array2D.Init(iWidth, iHeight)==false)
			return false;
	}

	iOffsetEndRow = iBytesPerLine - iWidth;
	pFloat = array2D.GetBuffer();
	pBits = pElements;

	for (y=0;y<iHeight;y++)
	{
		for (x=0;x<iWidth;x++)
		{
			*pFloat = (float)(*pBits)/255.0f;
			pBits++;
			pFloat++;
		}
		pBits += iOffsetEndRow;
	}
	return true;
}

bool CImage2D::ConvertToArray2DTripletFloatRGB(CArray2D<CTriplet<float> > &array2D) const
{
	int x, y, iOffsetEndRow, iBytesPerPixel;
	CImage2DPixel *pBits;
	CTriplet<float> *pTripletFloat;

	#ifdef IMAGE2D_RUNTIME_CHECK
	if (iBitsPerPixel!=24 && iBitsPerPixel!=32)
	{
		cerr<<"ERROR in CImage2D::ConvertToArray2DTripletFloatRGB(...): image is not RGB"<<endl;
		return false;
	}
	#endif

	if (array2D.GetSize()!=GetSize())
	{
		if (array2D.Init(iWidth, iHeight)==false)
			return false;
	}

	pTripletFloat = array2D.GetBuffer();
	iBytesPerPixel = iBitsPerPixel/8;
	iOffsetEndRow = iBytesPerLine - iWidth*iBytesPerPixel;
	pBits = pElements;

	for (y=0;y<iHeight;y++)
	{
		for (x=0;x<iWidth;x++)
		{
			pTripletFloat->Set((float)pBits[2], (float)pBits[1], (float)pBits[0]);
			*pTripletFloat /= 255.0f;
			pBits+=iBytesPerPixel;
			pTripletFloat++;
		}
		pBits += iOffsetEndRow;
	}
	return true;
}

bool CImage2D::ConvertToArray2DTripletFloatYUV(CArray2D<CTriplet<float> > &array2D) const
{
	int x, y, iOffsetEndRow, iBytesPerPixel;
	CImage2DPixel *pBits;
	CTriplet<float> *pTripletFloat;
	float fRed, fGreen, fBlue;

	#ifdef IMAGE2D_RUNTIME_CHECK
	if (iBitsPerPixel!=24 && iBitsPerPixel!=32)
	{
		cerr<<"ERROR in CImage2D::ConvertToArray2DTripletFloatYUV(...): image is not RGB"<<endl;
		return false;
	}
	#endif

	if (array2D.GetSize()!=GetSize())
	{
		if (array2D.Init(iWidth, iHeight)==false)
			return false;
	}

	pTripletFloat = array2D.GetBuffer();
	iBytesPerPixel = iBitsPerPixel/8;
	iOffsetEndRow = iBytesPerLine - iWidth*iBytesPerPixel;
	pBits = pElements;

	for (y=0;y<iHeight;y++)
	{
		for (x=0;x<iWidth;x++)
		{
			fRed = (float)pBits[2];
			fGreen = (float)pBits[1];
			fBlue = (float)pBits[0];

			pTripletFloat->x =  (0.299f*fRed + 0.587f*fGreen + 0.114f*fBlue)/255.0f; // Luminance Y
			pTripletFloat->y = (-0.147f*fRed - 0.289f*fGreen + 0.436f*fBlue)/255.0f; // Chrominance U
			pTripletFloat->z =  (0.615f*fRed - 0.515f*fGreen - 0.100f*fBlue)/255.0f; // Chrominance V

			pBits+=iBytesPerPixel;
			pTripletFloat++;
		}
		pBits += iOffsetEndRow;
	}
	return true;
}

bool CImage2D::ConvertToArray2DTripletFloatLAB(CArray2D<CTriplet<float> > &array2D) const
{
	int x, y, iOffsetEndRow, iBytesPerPixel;
	CImage2DPixel *pBits;
	CTriplet<float> *pTripletFloat;
	CTriplet<float> fReferenceWhiteXYZ, fPixelXYZ, fNormalizedXYZ, fTransformedXYZ;
	CImage2DFloatRGBPixel fRGBCorrected;
	float fRed, fGreen, fBlue;

	static float fGamma = 2.2f; // Gamma value for RGB correction (sRGB)
	static float fEpsilon = 0.008856f, fKappa = 903.3f; // Constants given by the CIE standards (for XYZ->Lab conversion)

	#ifdef IMAGE2D_RUNTIME_CHECK
	if (iBitsPerPixel!=24 && iBitsPerPixel!=32)
	{
		cerr<<"ERROR in CImage2D::ConvertToArray2DTripletFloatLAB(...): image is not RGB"<<endl;
		return false;
	}
	#endif

	if (array2D.GetSize()!=GetSize())
	{
		if (array2D.Init(iWidth, iHeight)==false)
			return false;
	}

	pTripletFloat = array2D.GetBuffer();
	iBytesPerPixel = iBitsPerPixel/8;
	iOffsetEndRow = iBytesPerLine - iWidth*iBytesPerPixel;
	pBits = pElements;

	// Get reference white in XYZ color space
	fRGBCorrected.Set(1.0f, 1.0f, 1.0f);

	// Convert to XYZ
	fReferenceWhiteXYZ.x = 0.576700f*fRGBCorrected.fRed + 0.297361f*fRGBCorrected.fGreen + 0.0270328f*fRGBCorrected.fBlue;
	fReferenceWhiteXYZ.y = 0.185556f*fRGBCorrected.fRed + 0.627355f*fRGBCorrected.fGreen + 0.0706879f*fRGBCorrected.fBlue;
	fReferenceWhiteXYZ.z = 0.188212f*fRGBCorrected.fRed + 0.075285f*fRGBCorrected.fGreen + 0.9912480f*fRGBCorrected.fBlue;

	for (y=0;y<iHeight;y++)
	{
		for (x=0;x<iWidth;x++)
		{
			fRed = (float)pBits[2];
			fGreen = (float)pBits[1];
			fBlue = (float)pBits[0];

			// The conversion methods were taken from
			// http://www.brucelindbloom.com/

			// RGB correction with gamma function
			fRGBCorrected.fRed = pow(fRed/255.0f, fGamma);
			fRGBCorrected.fGreen = pow(fGreen/255.0f, fGamma);
			fRGBCorrected.fBlue = pow(fBlue/255.0f, fGamma);

			// Convert to XYZ
			fPixelXYZ.x = 0.576700f*fRGBCorrected.fRed + 0.297361f*fRGBCorrected.fGreen + 0.0270328f*fRGBCorrected.fBlue;
			fPixelXYZ.y = 0.185556f*fRGBCorrected.fRed + 0.627355f*fRGBCorrected.fGreen + 0.0706879f*fRGBCorrected.fBlue;
			fPixelXYZ.z = 0.188212f*fRGBCorrected.fRed + 0.075285f*fRGBCorrected.fGreen + 0.9912480f*fRGBCorrected.fBlue;

			// Normalize with respect to reference white
			fNormalizedXYZ.x = fPixelXYZ.x/fReferenceWhiteXYZ.x;
			fNormalizedXYZ.y = fPixelXYZ.y/fReferenceWhiteXYZ.y;
			fNormalizedXYZ.z = fPixelXYZ.z/fReferenceWhiteXYZ.z;

			// Transform
			if (fNormalizedXYZ.x<=fEpsilon)
				fTransformedXYZ.x = (fKappa*fNormalizedXYZ.x+16.0f)/116.0f;
			else
				fTransformedXYZ.x = pow(fNormalizedXYZ.x, 1.0f/3.0f);
			if (fNormalizedXYZ.y<=fEpsilon)
				fTransformedXYZ.y = (fKappa*fNormalizedXYZ.y+16.0f)/116.0f;
			else
				fTransformedXYZ.y = pow(fNormalizedXYZ.y, 1.0f/3.0f);
			if (fNormalizedXYZ.z<=fEpsilon)
				fTransformedXYZ.z = (fKappa*fNormalizedXYZ.z+16.0f)/116.0f;
			else
				fTransformedXYZ.z = pow(fNormalizedXYZ.z, 1.0f/3.0f);

			// These is the standard scaling to obtain
			// L in [0,100], a in [-100,100] and b in [-100,100]
			// L = 116.0f*fTransformedXYZ.y-16.0f;
			// a = 500.0f*(fTransformedXYZ.x-fTransformedXYZ.y);
			// b = 200.0f*(fTransformedXYZ.y-fTransformedXYZ.z);

			// We actually apply the following scaling to have approximately
			// L in [0,1], a in [-0.5,0.5] and b in [-0.5,0.5]
			pTripletFloat->x = fTransformedXYZ.y; // L
			pTripletFloat->y = 2.5f*(fTransformedXYZ.x-fTransformedXYZ.y); // a
			pTripletFloat->z = fTransformedXYZ.y-fTransformedXYZ.z; // b

			pBits+=iBytesPerPixel;
			pTripletFloat++;
		}
		pBits += iOffsetEndRow;
	}

	return true;
}

bool CImage2D::ConvertToArray2DFloatRGBPixel(CArray2D<CImage2DFloatRGBPixel> &array2D) const
{
	int x, y, iOffsetEndRow, iBytesPerPixel;
	CImage2DPixel *pBits;
	CImage2DFloatRGBPixel *pPixel;

	#ifdef IMAGE2D_RUNTIME_CHECK
	if (iBitsPerPixel!=24 && iBitsPerPixel!=32)
	{
		cerr<<"ERROR in CImage2D::ConvertToArray2DFloatRGBPixel(...): image is not RGB"<<endl;
		return false;
	}
	#endif

	if (array2D.GetSize()!=GetSize())
	{
		if (array2D.Init(iWidth, iHeight)==false)
			return false;
	}

	pPixel = array2D.GetBuffer();
	iBytesPerPixel = iBitsPerPixel/8;
	iOffsetEndRow = iBytesPerLine - iWidth*iBytesPerPixel;
	pBits = pElements;

	for (y=0;y<iHeight;y++)
	{
		for (x=0;x<iWidth;x++)
		{
			pPixel->fRed = (float)pBits[2];
			pPixel->fGreen = (float)pBits[1];
			pPixel->fBlue = (float)pBits[0];
			*pPixel /= 255.0f;
			pBits+=iBytesPerPixel;
			pPixel++;
		}
		pBits += iOffsetEndRow;
	}
	return true;
}

void CImage2D::DrawLine(const CCouple<int> &p1, const CCouple<int> &p2, CImage2DPixel col)
{
	int x, y, xmin, xmax, ymin, ymax;
	float a, b;
	CImage2DPixel *pBits;

	#ifdef IMAGE2D_RUNTIME_CHECK
	if (iBitsPerPixel!=8)
	{
		cerr<<"ERROR in CImage2D::DrawLine(...): image is not 8-bit"<<endl;
		return;
	}
	#endif

	pBits = pElements;

	if (p1==p2)
	{
		if (p1.x>=0 && p1.x<iWidth && p1.y>=0 && p1.y<iHeight)
			pBits[p1.y*iBytesPerLine+p1.x] = col;
	}
	else if (p1.x!=p2.x)
	{
		a = (float)(p1.y-p2.y)/(p1.x-p2.x);
		b = (float)p1.y-a*p1.x;

		if (fabs(a)<=1.0f)
		{
			// Slope is greater in x than in y
			// Loop on x and compute corresponding y
			if (p1.x<p2.x) {xmin = p1.x; xmax = p2.x;}
			else {xmin = p2.x; xmax = p1.x;}

			xmin = max(xmin, 0);
			xmax = min(xmax, iWidth-1);

			for (x=xmin;x<=xmax;x++)
			{
				y = (int)(a*x+b);
				if (y>=0 && y<iHeight)
					pBits[y*iBytesPerLine+x]=col;
			}
		}
		else {
			// Slope is greater in y than in x
			// Loop on y and compute corresponding x
			if (p1.y<p2.y) {ymin = p1.y; ymax = p2.y;}
			else {ymin = p2.y; ymax = p1.y;}

			ymin = max(ymin, 0);
			ymax = min(ymax, iHeight-1);

			for (y=ymin;y<=ymax;y++)
			{
				x = (int)(((float)y-b)/a);
				if (x>=0 && x<iWidth)
					pBits[y*iBytesPerLine+x]=col;
			}
		}
	}
	else if (p1.x>=0 && p1.x<iWidth)
	{
		if (p1.y<p2.y) {ymin = p1.y; ymax = p2.y;}
		else {ymin = p2.y; ymax = p1.y;}

		ymin = max(ymin,0);
		ymax = min(ymax, iHeight-1);

		for (y=ymin;y<=ymax;y++)
			pBits[y*iBytesPerLine+p1.x]=col;
	}
}

void CImage2D::DrawFilledRectangle(const CCouple<int> &p1, const CCouple<int> &p2, CImage2DPixel col)
{
	CCouple<int> pmin, pmax;
	CImage2DPixel *pBits;
	int y, iRectangleWidth;

	#ifdef IMAGE2D_RUNTIME_CHECK
	if (iBitsPerPixel!=8)
	{
		cerr<<"ERROR in CImage2D::DrawFilledRectangle(...): image is not 8-bit"<<endl;
		return;
	}
	#endif

	pmin = coupleMax(coupleMin(p1, p2), CCouple<int>(0,0));
	pmax = coupleMin(coupleMax(p1, p2), CCouple<int>(iWidth-1, iHeight-1));
	iRectangleWidth = pmax.x - pmin.x + 1;

	pBits = pElements + pmin.y*iBytesPerLine + pmin.x;
	for (y=pmin.y;y<=pmax.y;y++)
	{
		memset(pBits, col, iRectangleWidth);
		pBits += iBytesPerLine;
	}
}

void CImage2D::DrawFilledCircle(const CCouple<int> &centre, int rayon, CImage2DPixel col)
{
	int y, y2, xmin, xmax, ymin, ymax;
	CImage2DPixel *pBits;

	#ifdef IMAGE2D_RUNTIME_CHECK
	if (iBitsPerPixel!=8)
	{
		cerr<<"ERROR in CImage2D::DrawFilledCircle(...): image is not 8-bit"<<endl;
		return;
	}
	#endif

	pBits = pElements;

	if (centre.y<rayon)
		ymin = -centre.y;
	else ymin = -rayon;

	if (centre.y+rayon>=iHeight)
		ymax = iHeight-centre.y-1;
	else ymax = rayon;

	for (y=ymin; y<=ymax; y++)
	{
		y2 = (int)sqrt((float)(rayon*rayon-y*y));
		if (centre.x<y2)
			xmin = -centre.x;
		else xmin = -y2;

		if (centre.x+y2>=iWidth)
			xmax = iWidth-centre.x-1;
		else xmax = y2;

		memset(pBits+(centre.y+y)*iBytesPerLine+(centre.x+xmin), col, xmax-xmin+1);
	}
}

void CImage2D::DrawLineRGB(const CCouple<int> &p1, const CCouple<int> &p2, const CImage2DByteRGBPixel &rgb)
{
	int x, y, xmin, xmax, ymin, ymax;
	float a, b;
	int iBytesPerPixel;
	CImage2DPixel *pBits;

	#ifdef IMAGE2D_RUNTIME_CHECK
	if (iBitsPerPixel!=24 && iBitsPerPixel!=32)
	{
		cerr<<"ERROR in CImage2D::DrawLineRGB(...): image is not RGB"<<endl;
		return;
	}
	#endif

	pBits = pElements;
	iBytesPerPixel = iBitsPerPixel/8;

	if (p1==p2)
	{
		if (p1.x>=0 && p1.x<iWidth && p1.y>=0 && p1.y<iHeight)
			*((CImage2DByteRGBPixel *)(pBits + p1.y*iBytesPerLine + p1.x*iBytesPerPixel)) = rgb;
	}
	else if (p1.x!=p2.x)
	{
		a = (float)(p1.y-p2.y)/(p1.x-p2.x);
		b = (float)p1.y-a*p1.x;

		if (fabs(a)<=1.0f)
		{
			// Slope is greater in x than in y
			// Loop on x and compute corresponding y
			if (p1.x<p2.x) {xmin = p1.x; xmax = p2.x;}
			else {xmin = p2.x; xmax = p1.x;}

			xmin = max(xmin, 0);
			xmax = min(xmax, iWidth-1);

			for (x=xmin;x<=xmax;x++)
			{
				y = (int)(a*x+b);
				if (y>=0 && y<iHeight)
					*((CImage2DByteRGBPixel *)(pBits + y*iBytesPerLine + x*iBytesPerPixel)) = rgb;
			}
		}
		else {
			// Slope is greater in y than in x
			// Loop on y and compute corresponding x
			if (p1.y<p2.y) {ymin = p1.y; ymax = p2.y;}
			else {ymin = p2.y; ymax = p1.y;}

			ymin = max(ymin, 0);
			ymax = min(ymax, iHeight-1);

			for (y=ymin;y<=ymax;y++)
			{
				x = (int)(((float)y-b)/a);
				if (x>=0 && x<iWidth)
					*((CImage2DByteRGBPixel *)(pBits + y*iBytesPerLine + x*iBytesPerPixel)) = rgb;
			}
		}
	}
	else if (p1.x>=0 && p1.x<iWidth)
	{
		if (p1.y<p2.y) {ymin = p1.y; ymax = p2.y;}
		else {ymin = p2.y; ymax = p1.y;}

		ymin = max(ymin,0);
		ymax = min(ymax, iHeight-1);

		for (y=ymin;y<=ymax;y++)
			*((CImage2DByteRGBPixel *)(pBits + y*iBytesPerLine + p1.x*iBytesPerPixel)) = rgb;
	}
}

void CImage2D::DrawFilledRectangleRGB(const CCouple<int> &p1, const CCouple<int> &p2, const CImage2DByteRGBPixel &rgb)
{
	CCouple<int> pmin, pmax;
	CImage2DPixel *pBits;
	int x, y, iBytesPerPixel, iOffsetEndRow;

	#ifdef IMAGE2D_RUNTIME_CHECK
	if (iBitsPerPixel!=24 && iBitsPerPixel!=32)
	{
		cerr<<"ERROR in CImage2D::DrawFilledRectangleRGB(...): image is not RGB"<<endl;
		return;
	}
	#endif

	pmin = coupleMax(coupleMin(p1, p2), CCouple<int>(0,0));
	pmax = coupleMin(coupleMax(p1, p2), CCouple<int>(iWidth-1, iHeight-1));

	iBytesPerPixel = iBitsPerPixel/8;
	iOffsetEndRow = iBytesPerLine - iWidth*iBytesPerPixel + (iWidth-1-pmax.x)*iBytesPerPixel + pmin.x*iBytesPerPixel;

	pBits = pElements + pmin.y*iBytesPerLine + pmin.x*iBytesPerPixel;
	for (y=pmin.y;y<=pmax.y;y++)
	{
		for (x=pmin.x; x<=pmax.x; x++)
		{
			*((CImage2DByteRGBPixel *)pBits) = rgb;
			pBits += iBytesPerPixel;
		}
		pBits += iOffsetEndRow;
	}
}

void CImage2D::DrawFilledCircleRGB(const CCouple<int> &centre, int rayon, const CImage2DByteRGBPixel &rgb)
{
	int x, y, y2, xmin, xmax, ymin, ymax, iBytesPerPixel;
	CImage2DPixel *pBits;

	#ifdef IMAGE2D_RUNTIME_CHECK
	if (iBitsPerPixel!=24 && iBitsPerPixel!=32)
	{
		cerr<<"ERROR in CImage2D::DrawFilledCircleRGB(...): image is not RGB"<<endl;
		return;
	}
	#endif

	iBytesPerPixel = iBitsPerPixel/8;

	if (centre.y<rayon)
		ymin = -centre.y;
	else ymin = -rayon;

	if (centre.y+rayon>=iHeight)
		ymax = iHeight-centre.y-1;
	else ymax = rayon;

	for (y=ymin; y<=ymax; y++)
	{
		y2 = (int)sqrt((float)(rayon*rayon-y*y));
		if (centre.x<y2)
			xmin = -centre.x;
		else xmin = -y2;

		if (centre.x+y2>=iWidth)
			xmax = iWidth-centre.x-1;
		else xmax = y2;

		pBits = pElements+(centre.y+y)*iBytesPerLine+(centre.x+xmin)*iBytesPerPixel;
		for (x=xmin; x<=xmax; x++)
		{
			*((CImage2DByteRGBPixel *)pBits) = rgb;
			pBits+=iBytesPerPixel;
		}
	}
}

void CImage2D::Clear(CImage2DPixel intensity)
{
	CImage2DPixel *pBits;
	int y;

	#ifdef IMAGE2D_RUNTIME_CHECK
	if (iBitsPerPixel!=8)
	{
		cerr<<"ERROR in CImage2D::Clear(...): image is not 8-bit"<<endl;
		return;
	}
	#endif

	pBits = pElements;
	for (y=0;y<iHeight;y++)
	{
		memset(pBits, intensity, iWidth);
		pBits+=iBytesPerLine;
	}
}

void CImage2D::ClearRGB(const CImage2DByteRGBPixel &rgb)
{
	#ifdef IMAGE2D_RUNTIME_CHECK
	if (iBitsPerPixel!=24 && iBitsPerPixel!=32)
	{
		cerr<<"ERROR in CImage2D::ClearRGB(...): image is not RGB"<<endl;
		return;
	}
	#endif

	if (iBitsPerPixel==32)
	{
		RGBQUAD rgbQuad;
		RGBQUAD *pPixels = (RGBQUAD *)pElements;
		int i, iNbPixels = iWidth*iHeight;

		rgbQuad.rgbBlue = rgb.byteBlue;
		rgbQuad.rgbGreen = rgb.byteGreen;
		rgbQuad.rgbRed = rgb.byteRed;
		rgbQuad.rgbReserved = 0;

		for (i=0;i<iNbPixels;i++)
		{
			*pPixels = rgbQuad;
			pPixels++;
		}
	}
	else if (iBitsPerPixel==24)
	{
		int x, y, iOffsetEndRow;
		iOffsetEndRow = iBytesPerLine - 3*iWidth;
		CImage2DByteRGBPixel *pPixels = (CImage2DByteRGBPixel *)pElements;

		for (y=0;y<iHeight;y++)
		{
			for (x=0;x<iWidth;x++)
			{
				*pPixels = rgb;
				pPixels++;
			}
			pPixels = (CImage2DByteRGBPixel *)(((unsigned char *)pPixels)+iOffsetEndRow);
		}
	}
}

void CImage2D::Flip()
{
	int y;
	CImage2DPixel *pBitsStart, *pBitsEnd, *pBitsTemp;

	pBitsStart = pElements;
	pBitsEnd = pElements + (iHeight-1)*iBytesPerLine;

	pBitsTemp = new CImage2DPixel[iBytesPerLine];
	for (y=0;y<iHeight/2;y++)
	{
		memcpy(pBitsTemp, pBitsStart, iBytesPerLine);
		memcpy(pBitsStart, pBitsEnd, iBytesPerLine);
		memcpy(pBitsEnd, pBitsTemp, iBytesPerLine);

		pBitsStart += iBytesPerLine;
		pBitsEnd -= iBytesPerLine;
	}
	delete[] pBitsTemp;
}

// BMP Encoder/Decoder
bool CImage2D::DecodeBMP(FILE *pFile)
{
	BITMAPFILEHEADER bmpFileHeader;
	BITMAPINFOHEADER bmpInfoHeader;
	CArray1D<RGBQUAD> arrayPalette;
	unsigned int off;
	bool bHasPalette = false;
	int iBitsPerPixel;
	size_t sizeBmpFileHeaderMin = (size_t)14;

	if (pFile == NULL)
		return false;

	off = ftell(pFile);

	if (fread(&bmpFileHeader,min(sizeBmpFileHeaderMin,sizeof(bmpFileHeader)),1,pFile)==0)
	{
		cerr<<"ERROR in CImage2D::DecodeBMP(...): cannot read BMP file header"<<endl;
		return false;
	}

	if (bmpFileHeader.bfType != 0x4D42) // (bmpFileHeader.bfType != BFT_BITMAP) // Do we have a RC HEADER ?
	{
        bmpFileHeader.bfOffBits = 0L;
        fseek(pFile,off,SEEK_SET);
    }

    if (fread(&bmpInfoHeader,sizeof(BITMAPINFOHEADER),1,pFile)==0)
	{
		cerr<<"ERROR in CImage2D::DecodeBMP(...): cannot read BMP info header"<<endl;
		return false;
	}

	if (bmpInfoHeader.biSize!=sizeof(BITMAPINFOHEADER))
	{
		cerr<<"ERROR in CImage2D::DecodeBMP(...): unusual size of BMP file header"<<endl;
		return false;
	}

	if (bmpInfoHeader.biBitCount==8)
	{
		// Compute size of supposed palette in bytes
		unsigned int uiNbPaletteBytes =
			bmpFileHeader.bfOffBits - sizeof(BITMAPINFOHEADER) - min(sizeBmpFileHeaderMin,sizeof(bmpFileHeader));

		if (uiNbPaletteBytes==0)
		{
			bHasPalette = false;
			iBitsPerPixel = 8;
		}
		else if (uiNbPaletteBytes%4==0)
		{
			unsigned int uiNbColorsInPalette;

			uiNbColorsInPalette = uiNbPaletteBytes/4;
			arrayPalette.Init(uiNbColorsInPalette);
			bHasPalette = true;
			iBitsPerPixel = 24;

            if (fread(arrayPalette.GetBuffer(), sizeof(RGBQUAD), uiNbColorsInPalette, pFile)!=uiNbColorsInPalette)
			{
				cerr<<"ERROR in CImage2D::DecodeBMP(...): cannot read all colors in palette. Unexpected end-of-file ?"<<endl;
				return false;
			}
		}
		else {
			cerr<<"ERROR in CImage2D::DecodeBMP(...): unusual number of bytes between info header and image data"<<endl;
			return false;
		}
	}
	else iBitsPerPixel = bmpInfoHeader.biBitCount;

	if (Create(bmpInfoHeader.biWidth, abs(bmpInfoHeader.biHeight), iBitsPerPixel)==false)
	{
		cerr<<"ERROR in CImage2D::DecodeBMP(...): cannot allocate memory"<<endl;
		return false;
	}

	//if (bmpFileHeader.bfOffBits != 0L)
	//	fseek(pFile, off + bmpFileHeader.bfOffBits, SEEK_SET);

	if (bmpInfoHeader.biCompression==BI_RGB)
	{
		if ((bmpInfoHeader.biBitCount==8 && bHasPalette==false) || bmpInfoHeader.biBitCount==24 || bmpInfoHeader.biBitCount==32)
		{
			unsigned int iSize = iBytesPerLine*iHeight;
			if (fread(pElements, 1, iSize, pFile)!=iSize)
			{
				cerr<<"ERROR in CImage2D::DecodeBMP(...): cannot read all pixels. Unexpected end-of-file ?"<<endl;
				return false;
			}
		}
		else if (bmpInfoHeader.biBitCount==8 && bHasPalette==true)
		{
			int x, y;
			CImage2DPixel *pPixels;
			CArray1D<unsigned char> rowIndices;
			int iOffsetEndRow;
			unsigned int iBytesPerRowSource;

			pPixels = pElements;
			iOffsetEndRow = iBytesPerLine - 3*iWidth;
			if ((iWidth+4)%4!=0)
				iBytesPerRowSource = iWidth + 4-(iWidth+4)%4;
			else
				iBytesPerRowSource = iWidth;

			rowIndices.Init(iBytesPerRowSource);

			for (y=0; y<iHeight; y++)
			{
				if (fread(rowIndices.GetBuffer(), 1, iBytesPerRowSource, pFile)!=iBytesPerRowSource)
				{
					cerr<<"ERROR in CImage2D::DecodeBMP(...): cannot read all pixels. Unexpected end-of-file ?"<<endl;
					return false;
				}
				for (x=0; x<iWidth; x++)
				{
					*(CImage2DByteRGBPixel *)pPixels = *((CImage2DByteRGBPixel *)(&(arrayPalette[rowIndices[x]])));
					pPixels += 3;
				}
				pPixels += iOffsetEndRow;
			}
		}
		else {
			cerr<<"ERROR in CImage2D::DecodeBMP(...): unsupported RGB depth"<<endl;
			return false;
		}
	}
	else if (bmpInfoHeader.biCompression==BI_RLE8)
	{
		unsigned char status_byte = 0;
		unsigned char second_byte = 0;
		int scanline = 0;
		int bits = 0;
		bool bContinue = true;
		CImage2DPixel *pRow;
		CArray1D<unsigned char> arrayBufferRead;
        int x;
        CImage2DPixel *pPixels;
        unsigned char *pReadBytes;

		const int RLE_COMMAND     = 0;
		const int RLE_ENDOFLINE   = 0;
		const int RLE_ENDOFBITMAP = 1;
		const int RLE_DELTA       = 2;

		while (bContinue==true)
		{
			if (fread(&status_byte, sizeof(unsigned char), 1, pFile)!=1)
			{
                cerr<<"ERROR in CImage2D::DecodeBMP(...): cannot read status byte in RLE encoding. Unexpected end-of-file ?"<<endl;
                return false;
			}
			switch (status_byte)
			{
				case RLE_COMMAND:
					if (fread(&status_byte, sizeof(unsigned char), 1, pFile)!=1)
                    {
                        cerr<<"ERROR in CImage2D::DecodeBMP(...): cannot read command byte in RLE encoding. Unexpected end-of-file ?"<<endl;
                        return false;
                    }
					switch (status_byte) {
						case RLE_ENDOFLINE:
							bits = 0;
							scanline++;
							if (scanline==iHeight)
								bContinue = false;
							break;
						case RLE_ENDOFBITMAP:
							bContinue = false;
							break;
						case RLE_DELTA:
						{
							// read the delta values
							unsigned char delta_x;
							unsigned char delta_y;
							if (fread(&delta_x, sizeof(unsigned char), 1, pFile)!=1)
                            {
                                cerr<<"ERROR in CImage2D::DecodeBMP(...): cannot read DeltaX byte in RLE encoding. Unexpected end-of-file ?"<<endl;
                                return false;
                            }
							if (fread(&delta_y, sizeof(unsigned char), 1, pFile)!=1)
                            {
                                cerr<<"ERROR in CImage2D::DecodeBMP(...): cannot read DeltaX byte in RLE encoding. Unexpected end-of-file ?"<<endl;
                                return false;
                            }
							// apply them
							bits     += delta_x;
							scanline += delta_y;
							break;
						}
						default:
							pRow = pElements + scanline*iBytesPerLine;
							arrayBufferRead.Init(status_byte);

							if (fread(arrayBufferRead.GetBuffer(), sizeof(unsigned char), status_byte, pFile)!=status_byte)
                            {
                                cerr<<"ERROR in CImage2D::DecodeBMP(...): cannot read "
                                    <<(unsigned int)status_byte<<" required bytes. Unexpected end-of-file ?"<<endl;
                                return false;
                            }

                            // If image has a palette, get RGB colors from indices
                            // Otherwise, just copy read bytes to the image data
                            if (iBitsPerPixel==24 && bHasPalette==true)
                            {
                                pPixels = pRow + bits;
                                pReadBytes = arrayBufferRead.GetBuffer();
                                for (x=0; x<status_byte; x++)
                                {
                                    *(CImage2DByteRGBPixel *)pPixels = *((CImage2DByteRGBPixel *)(&(arrayPalette[*pReadBytes])));
                                    pReadBytes++;
                                    pPixels += 3;
                                }
                                bits += status_byte*3;
                            }
                            else {
                                memcpy(pRow + bits, arrayBufferRead.GetBuffer(), status_byte);
                                bits += status_byte;
                            }
							// align run length to even number of bytes
							if ((status_byte & 1) == 1)
								if (fread(&second_byte, sizeof(unsigned char), 1, pFile)!=1)
                                {
                                    cerr<<"ERROR in CImage2D::DecodeBMP(...): cannot read second byte. Unexpected end-of-file ?"<<endl;
                                    return false;
                                }
							// bits += status_byte;
							break;
					};
					break;
				default:
					pRow = pElements + scanline*iBytesPerLine;
					if (fread(&second_byte, sizeof(unsigned char), 1, pFile)!=1)
                    {
                        cerr<<"ERROR in CImage2D::DecodeBMP(...): cannot read second byte. Unexpected end-of-file ?"<<endl;
                        return false;
                    }

                    if (bHasPalette==true)
                    {
                        for (int i = 0; i < status_byte; i++)
                        {
                            if (bits<iBytesPerLine)
                            {
                                *((CImage2DByteRGBPixel *)(pRow + bits)) = *((CImage2DByteRGBPixel *)(&(arrayPalette[second_byte])));
                                // *(pRow + bits) = second_byte;
                                bits+=3;
                            }
                            else
                                bContinue = false;
                        }

                    }
                    else {
                        for (int i = 0; i < status_byte; i++)
                        {
                            if (bits<iBytesPerLine)
                            {
                                *(pRow + bits) = second_byte;
                                bits++;
                            }
                            else
                                bContinue = false;
                        }
                    }
					break;
			};
		}
	}
	else {
		cerr<<"ERROR in CImage2D::DecodeBMP(...): unsupported compression mode"<<endl;
		return false;
	}

	// If biHeight fields is positive, file is a a bottom-up bitmap, so image should be flipped
	if (bmpInfoHeader.biHeight>0)
		Flip();

	return true;
}

bool CImage2D::EncodeBMP(FILE *pFile) const
{
	BITMAPFILEHEADER bmpFileHeader;
	BITMAPINFOHEADER bmpInfoHeader;
	unsigned int iSizeTotal;
	unsigned int iSizePalette;

	if (pFile==NULL)
		return false;

	if (iBitsPerPixel==8)
		iSizePalette = sizeof(RGBQUAD)*256;
	else if (iBitsPerPixel==24 || iBitsPerPixel==32)
		iSizePalette = 0;
	else {
		cerr<<"ERROR in CImage2D::EncodeBMP(...): Unsupported number of bits per pixel"<<endl;
		return false;
	}
	iSizeTotal = sizeof(BITMAPFILEHEADER)+sizeof(BITMAPINFOHEADER)+iSizePalette+iBytesPerLine*iHeight;

	// Set file header fields
	bmpFileHeader.bfType = 0x4D42; // 'BM'
	bmpFileHeader.bfSize = iSizeTotal;
	bmpFileHeader.bfReserved1 = 0;
	bmpFileHeader.bfReserved2 = 0;
	bmpFileHeader.bfOffBits = sizeof(BITMAPFILEHEADER)+sizeof(BITMAPINFOHEADER)+iSizePalette;

	// Set info header fields
    bmpInfoHeader.biSize = sizeof(BITMAPINFOHEADER);
    bmpInfoHeader.biWidth = iWidth;
    bmpInfoHeader.biHeight = -iHeight;
    bmpInfoHeader.biPlanes = 1;
    bmpInfoHeader.biBitCount = (unsigned short int)iBitsPerPixel;
    bmpInfoHeader.biCompression = BI_RGB;
    bmpInfoHeader.biSizeImage = iBytesPerLine*iHeight;
    bmpInfoHeader.biXPelsPerMeter = 96;
    bmpInfoHeader.biYPelsPerMeter = 96;
    if (iBitsPerPixel==8)
		bmpInfoHeader.biClrUsed = 256;
	else
		bmpInfoHeader.biClrUsed = 0;
    bmpInfoHeader.biClrImportant = 0;

	// Write file and info headers
	fwrite(&bmpFileHeader, sizeof(BITMAPFILEHEADER), 1, pFile);
	fwrite(&bmpInfoHeader, sizeof(BITMAPINFOHEADER), 1, pFile);

	// Write palette if image is 8-bit indexed
	if (iSizePalette!=0)
	{
		RGBQUAD *pPalette;
		pPalette = (RGBQUAD *)(CImage2D::vectBitsBitmapInfo.GetBuffer() + sizeof(BITMAPINFOHEADER));
		fwrite(pPalette, iSizePalette, 1, pFile);
	}

	// Write pixel data
	fwrite(pElements, 1, iBytesPerLine*iHeight, pFile);

	return true;
}

// JPEG Encoder/Decoder
#ifdef IMAGE2D_SUPPORT_JPG

struct jpg_error_mgr {
	struct jpeg_error_mgr pub;	// "public" fields
	jmp_buf setjmp_buffer;		// for return to caller
	char* buffer;				// error message
};
typedef jpg_error_mgr *jpg_error_ptr;

// Here's the routine that will replace the standard error_exit method:
void ima_jpeg_error_exit(j_common_ptr cinfo)
{
	// cinfo->err really points to a my_error_mgr struct, so coerce pointer
	jpg_error_ptr myerr = (jpg_error_ptr) cinfo->err;

	// Create the message
	myerr->pub.format_message(cinfo, myerr->buffer);

	// Send it to stderr, adding a newline
	// Return control to the setjmp point
	longjmp(myerr->setjmp_buffer, 1);
}

bool CImage2D::DecodeJPG(FILE *pFile)
{
	char strErrorBuffer[2000];
	CImage2DPixel *pRow;

	// This structure contains the JPEG decompression parameters and pointers to working space
	struct jpeg_decompress_struct cinfo;

	// We use our private extension JPEG error handler
	struct jpg_error_mgr jerr;

	// Output row buffer
	JSAMPARRAY buffer;

	if (pFile==NULL)
	{
		cerr<<"ERROR in CImage2D::DecodeJPG(...): file pointer parameter is NULL"<<endl;
		return false;
	}

	// Set error buffer
	jerr.buffer = strErrorBuffer;

	// Set up the normal JPEG error routines, then override error_exit.
	cinfo.err = jpeg_std_error(&jerr.pub);
	jerr.pub.error_exit = ima_jpeg_error_exit;

	// Establish the setjmp return context for my_error_exit to use.
	if (setjmp(jerr.setjmp_buffer))
	{
		// If we get here, the JPEG code has signaled an error.
		// We need to clean up the JPEG object, close the input file, and return.
		jpeg_destroy_decompress(&cinfo);
		return false;
	}

	// Initialize decompression structure
	jpeg_create_decompress(&cinfo);

	// Specify data source
	jpeg_stdio_src(&cinfo, pFile);

	// Read file parameters
	jpeg_read_header(&cinfo, TRUE);

	// Set the scale
	// cinfo.scale_denom = 1;

	// Start decompressor
	jpeg_start_decompress(&cinfo);

	// Create the image using output dimensions
	if (Create(cinfo.output_width, cinfo.output_height, 8*cinfo.output_components)==false)
		longjmp(jerr.setjmp_buffer, 1); // Check if the image has been created

	// Make a one-row-high sample array that will go away when done with image
	buffer = (*cinfo.mem->alloc_sarray)((j_common_ptr) &cinfo, JPOOL_IMAGE, iBytesPerLine, 1);

	// Read all scanlines
	pRow = pElements;
	while (cinfo.output_scanline < cinfo.output_height)
	{
		jpeg_read_scanlines(&cinfo, buffer, 1);

		// CMYK -> RGB
		if ((cinfo.num_components==4)&&(cinfo.quantize_colors==FALSE))
		{
			unsigned char k,*dst,*src;

			dst = pRow;
			src = buffer[0];

			for (long x3=0,x4=0; x3<iBytesPerLine && x4<iBytesPerLine; x3+=3, x4+=4)
			{
				k = src[x4+3];
				dst[x3]   = (unsigned char)((k * src[x4+2])/255);
				dst[x3+1] = (unsigned char)((k * src[x4+1])/255);
				dst[x3+2] = (unsigned char)((k * src[x4+0])/255);
			}
		}
		else
			memcpy(pRow, buffer[0], iBytesPerLine);

		pRow += iBytesPerLine;
	}

	// Finish decompression
	jpeg_finish_decompress(&cinfo);

	// Swap red and blue components
	// Not necessary if swapped red and blue definition in jmorecfg.h
	if ((cinfo.num_components==3)&&(cinfo.quantize_colors==FALSE))
	{
		CImage2DPixel *pBits;
		CImage2DPixel byteTemp;
		int x, y, iOffsetEndRow, iBytesPerPixel;

		pBits = pElements;
		iBytesPerPixel = 3;
		iOffsetEndRow = iBytesPerLine - iWidth*iBytesPerPixel;
		for (y=0; y<iHeight; y++)
		{
			for (x=0; x<iWidth; x++)
			{
				byteTemp = pBits[2];
				pBits[2] = pBits[0];
				pBits[0] = byteTemp;
				pBits += iBytesPerPixel;
			}
			pBits += iOffsetEndRow;
		}
	}

	// Finish decompression ?
	// jpeg_finish_decompress(&cinfo); // Error generated if called

	// Release JPEG decompression structure
	jpeg_destroy_decompress(&cinfo);

	return true;
}

bool CImage2D::EncodeJPG(FILE *pFile) const
{
	char strErrorBuffer[2000];
	CImage2DPixel *pRow;

	// This structure contains the JPEG compression parameters and pointers to working space
	struct jpeg_compress_struct cinfo;

	// We use our private extension JPEG error handler
	struct jpg_error_mgr jerr;

	// Output row buffer
	JSAMPARRAY buffer;

	if (pFile==NULL)
	{
		cerr<<"ERROR in CImage2D::EncodeJPG(...): file pointer parameter is NULL"<<endl;
		return false;
	}

	if (iBitsPerPixel!=8 && iBitsPerPixel!=24)
	{
		cerr<<"ERROR in CImage2D::EncodeJPG(...): only 8 or 24 bits images can be saved in JPEG format"<<endl;
		return false;
	}

	// Set error buffer
	jerr.buffer = strErrorBuffer;

	// Set up the normal JPEG error routines, then override error_exit
	cinfo.err = jpeg_std_error(&jerr.pub);
	jerr.pub.error_exit = ima_jpeg_error_exit;

	// Establish the setjmp return context for my_error_exit to use
	if (setjmp(jerr.setjmp_buffer))
	{
		// If we get here, the JPEG code has signaled an error.
		// Clean up the JPEG object, close the input file, and return.
		strcpy(strErrorBuffer, jerr.buffer);
		jpeg_destroy_compress(&cinfo);
		return false;
	}

	// Initialize compression structure
	jpeg_create_compress(&cinfo);

	// Specify data destination
	jpeg_stdio_dest(&cinfo, pFile);

	// Set basic properties of image (size and number of components per pixel)
	cinfo.image_width = iWidth;
	cinfo.image_height = iHeight;

	// Set number of components per pixel, whether the image is 8-bit grayscale or 24-bit RGB
	if (iBitsPerPixel==8)
	{
		cinfo.input_components = 1;
		cinfo.in_color_space = JCS_GRAYSCALE;
	}
	else {
		cinfo.input_components = 3;
		cinfo.in_color_space = JCS_RGB;
	}

	// Set default compression parameters.
	jpeg_set_defaults(&cinfo);

	// Set quality parameter
	jpeg_set_quality(&cinfo, 100, 0);

	cinfo.density_unit = 1;
	cinfo.X_density = 96;
	cinfo.Y_density = 96;

	// Start compressor
	// TRUE ensures that we will write a complete interchange-JPEG file.
	jpeg_start_compress(&cinfo, TRUE);

	// Make a one-row-high sample array that will go away when done with image
	// "8+iBytesPerLine" fix heap deallocation problem during debug ???
	buffer = (*cinfo.mem->alloc_sarray)
		((j_common_ptr) &cinfo, JPOOL_IMAGE, 8+iBytesPerLine, 1);

	// Write all scanlines
	pRow = pElements;
	while (cinfo.next_scanline < cinfo.image_height)
	{
		memcpy(buffer[0], pRow, iBytesPerLine);

		// Swap red and blue components
		// Not necessary if swapped red and blue definition in jmorecfg.h
		if (iBitsPerPixel==24)
		{
			CImage2DPixel *pBits;
			CImage2DPixel byteTemp;
			int x, iBytesPerPixel;

			pBits = (CImage2DPixel *)(buffer[0]);
			iBytesPerPixel = 3;

			for (x=0; x<iWidth; x++)
			{
				byteTemp = pBits[2];
				pBits[2] = pBits[0];
				pBits[0] = byteTemp;
				pBits += iBytesPerPixel;
			}
		}
		jpeg_write_scanlines(&cinfo, buffer, 1);
		pRow += iBytesPerLine;
	}

	// Finish compression
	jpeg_finish_compress(&cinfo);

	// Release JPEG compression structure
	jpeg_destroy_compress(&cinfo);

	return true;
}
#endif // IMAGE2D_SUPPORT_JPG

#ifdef IMAGE2D_SUPPORT_PNG

void user_read_data(png_structp png_ptr, png_bytep data, png_size_t length)
{
	FILE *pFile = (FILE *)png_get_io_ptr(png_ptr);
	if (fread(data, 1, length, pFile)!=length)
		png_error(png_ptr, "Read error");
}

void user_write_data(png_structp png_ptr, png_bytep data, png_size_t length)
{
	FILE *pFile = (FILE *)png_get_io_ptr(png_ptr);
	if (fwrite(data, 1, length, pFile)!=length)
		png_error(png_ptr, "Write error");
}

void user_flush_data(png_structp png_ptr)
{
	FILE *pFile = (FILE *)png_get_io_ptr(png_ptr);
	if (fflush(pFile)==0)
		png_error(png_ptr, "Flush error");
}

void user_error_fn(png_structp png_ptr, png_const_charp error_msg)
{
	strncpy((char*)png_get_error_ptr(png_ptr),error_msg,255);
	longjmp(png_jmpbuf(png_ptr), 1);
}

bool CImage2D::DecodePNG(FILE *pFile)
{
	png_structp png_ptr;
	png_infop info_ptr;
	png_uint_32 width, height;
	int bit_depth, color_type, interlace_type;
	png_bytep *row_pointers;
	int row;

	png_ptr = NULL;
	info_ptr = NULL;
	row_pointers = NULL;

	if (pFile==NULL)
	{
		cerr<<"ERROR in CImage2D::DecodePNG(...): file pointer parameter is NULL"<<endl;
		return false;
	}

	// Create and initialize the png_struct with the desired error handler functions
	png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
	if (png_ptr == NULL)
	{
		cerr<<"ERROR in CImage2D::DecodePNG(...): cannot create PNG structure"<<endl;
		return false;
	}

	// Allocate and initialize the image information data
	info_ptr = png_create_info_struct(png_ptr);
	if (info_ptr == NULL)
	{
		cerr<<"ERROR in CImage2D::DecodePNG(...): cannot create PNG info structure"<<endl;
		png_destroy_read_struct(&png_ptr, (png_infopp)NULL, (png_infopp)NULL);
		return false;
	}

	// Set error handling
	if (setjmp(png_jmpbuf(png_ptr)))
	{
		cerr<<"ERROR in CImage2D::DecodePNG(...): cannot set error handling"<<endl;

		// Free all of the memory associated with the png_ptr and info_ptr
		if (row_pointers!=NULL)
			delete[] row_pointers;
		png_destroy_read_struct(&png_ptr, &info_ptr, (png_infopp)NULL);

		return false;
	}

	// Set up the input control
	png_set_read_fn(png_ptr, pFile, (png_rw_ptr)user_read_data);

	// Read header info
	png_read_info(png_ptr, info_ptr);
	png_get_IHDR(png_ptr, info_ptr, &width, &height, &bit_depth, &color_type, &interlace_type, NULL, NULL);

	// Tell libpng to strip 16 bit/color files down to 8 bits/color
	png_set_strip_16(png_ptr) ;

	// Extract multiple pixels with bit depths of 1, 2, and 4 from a single
	// byte into separate bytes (useful for paletted and grayscale images).
	png_set_packing(png_ptr);

	// Scale grayscale values to the range 0..255
	if (color_type == PNG_COLOR_TYPE_GRAY)
		png_set_expand(png_ptr);

	if (color_type == PNG_COLOR_TYPE_PALETTE)
		png_set_palette_to_rgb(png_ptr);

	if (color_type == PNG_COLOR_TYPE_GRAY_ALPHA)
		png_set_gray_to_rgb(png_ptr);

	// Flip the RGB pixels to BGR (or RGBA to BGRA)
	if (color_type & PNG_COLOR_MASK_COLOR)
		png_set_bgr(png_ptr);

	png_read_update_info(png_ptr, info_ptr);

	png_get_IHDR(png_ptr, info_ptr, &width, &height, &bit_depth, &color_type, &interlace_type, NULL, NULL);

	Create(width, height, bit_depth*png_get_channels(png_ptr, info_ptr));

	// Create the array of pointers to image data
	row_pointers = (png_bytep *)new png_bytep[height];
	if (row_pointers == NULL)
	{
		cerr<<"ERROR in CImage2D::DecodePNG(...): cannot allocate row pointers array"<<endl;
		png_destroy_read_struct(&png_ptr, (png_infopp)NULL, (png_infopp)NULL);
		return false;
	}

	for (row = 0; row < (int)height; row++)
		row_pointers[row] = (png_byte *)pElements + row*iBytesPerLine;

	// Read the entire image in one go
	png_read_image(png_ptr, row_pointers);

	delete[] row_pointers;

	// Read the rest of the file, getting any additional chunks in info_ptr
	png_read_end(png_ptr, info_ptr);

	// Clean up after the read, and free any memory allocated
	png_destroy_read_struct(&png_ptr, &info_ptr, (png_infopp)NULL);

	return true;
}

bool CImage2D::EncodePNG(FILE *pFile) const
{
	png_struct *png_ptr;
	png_info *info_ptr;
	unsigned char *row_pointers = NULL;
	CImage2DPixel *pRow;
	int y, pass, number_passes;
	int bit_depth;
	int color_type;
	int interlace_type;
	int compression_type;
	int filter_method;

	if (pFile==NULL)
	{
		cerr<<"ERROR in CImage2D::EncodePNG(...): file pointer parameter is NULL"<<endl;
		return false;
	}

	// Create and initialize the png_struct with the desired error handler functions
	png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING,(void *)NULL,NULL,NULL);
	if (png_ptr == NULL)
	{
		cerr<<"ERROR in CImage2D::EncodePNG(...): cannot create PNG structure"<<endl;
		return false;
	}

	// Allocate and initialize the image information data
	info_ptr = png_create_info_struct(png_ptr);
	if (info_ptr == NULL)
	{
		png_destroy_write_struct(&png_ptr, (png_infopp)NULL);
		cerr<<"ERROR in CImage2D::EncodePNG(...): cannot create PNG info structure"<<endl;
		return false;
	}

	// Set error handling
	if (setjmp(png_jmpbuf(png_ptr)))
	{
		cerr<<"ERROR in CImage2D::EncodePNG(...): cannot set error handling"<<endl;

		// Free all of the memory associated with the png_ptr and info_ptr
		png_destroy_write_struct(&png_ptr, (png_infopp)&info_ptr);
		return false;
	}

	// Use custom I/O functions
	png_set_write_fn(png_ptr,pFile,(png_rw_ptr)user_write_data,(png_flush_ptr)user_flush_data);

	bit_depth = 8;

	// 0 = gray, 1 = indexed (with palette), 2 = RGB, 4 = RGBA
	if (iBitsPerPixel==8)
		color_type = PNG_COLOR_TYPE_GRAY;
	else if (iBitsPerPixel==24)
		color_type = PNG_COLOR_TYPE_RGB;
	else if (iBitsPerPixel==32)
		color_type = PNG_COLOR_TYPE_RGB_ALPHA;
    else {
		cerr<<"ERROR in CImage2D::EncodePNG(...): number of bits/pixel is not valid"<<endl;

		// Free all of the memory associated with the png_ptr and info_ptr
		png_destroy_write_struct(&png_ptr, (png_infopp)&info_ptr);
		return false;
    }

	compression_type = 0;
	filter_method = 0;
	interlace_type = PNG_INTERLACE_NONE;

	// Set the file information
	png_set_IHDR(png_ptr, info_ptr, iWidth, iHeight,
		bit_depth, color_type, interlace_type, compression_type, filter_method);

	// Set background
	png_color_16 image_background={ 0, 255, 255, 255, 0 };
	png_set_bKGD(png_ptr, info_ptr, &image_background);

	// Set metrics
	png_set_pHYs(png_ptr, info_ptr, 96, 96, PNG_RESOLUTION_METER);

	// Write the file information
	png_write_info(png_ptr, info_ptr);

	// We will write one row at a time
	row_pointers = new unsigned char[10+iBytesPerLine];

	// Interlace handling
	number_passes = png_set_interlace_handling(png_ptr);

	for (pass = 0; pass < number_passes; pass++)
	{
		// Write all rows
		pRow = pElements;
		for (y=0; y<iHeight; y++)
		{
			memcpy(row_pointers, pRow, iBytesPerLine);

			//HACK BY OP
			if (color_type == PNG_COLOR_TYPE_RGB)
			{
				CImage2DPixel *pBits;
				CImage2DPixel byteTemp;
				int x, iBytesPerPixel;

				pBits = (CImage2DPixel *)row_pointers;
				iBytesPerPixel = 3;

				for (x=0; x<iWidth; x++)
				{
					byteTemp = pBits[2];
					pBits[2] = pBits[0];
					pBits[0] = byteTemp;
					pBits += iBytesPerPixel;
				}
			}
			png_write_row(png_ptr, row_pointers);
			pRow += iBytesPerLine;
		}
	}

	delete [] row_pointers;

	// Finish writing rest of the file
	png_write_end(png_ptr, info_ptr);

	// Free palette if needed
	//if (info_ptr->palette!=NULL)
	//	delete[] info_ptr->palette;

	// Clean up after writing, and free any memory allocated
	png_destroy_write_struct(&png_ptr, (png_infopp)&info_ptr);

	return true;
}

#endif // IMAGE2D_SUPPORT_PNG
