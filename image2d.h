#ifndef _IMAGE2D_H
#define _IMAGE2D_H

#include <algorithm>

#include "arrayndfloat.h"
#include "couple.h"
#include "triplet.h"

#if (defined(_WIN32) || defined(WIN32))
#ifdef _MSC_VER
#define NOMINMAX
#endif
#include <windows.h>
#endif

// Supported image formats (in addition to BMP)
// Implies linking with libjpeg
#define IMAGE2D_SUPPORT_JPG

// Implies linking with libpng and zlib
#define IMAGE2D_SUPPORT_PNG

// #define IMAGE2D_RUNTIME_CHECK

typedef unsigned char CImage2DPixel;

class CImage2DByteRGBPixel
{
  public:
	unsigned char byteBlue;
	unsigned char byteGreen;
	unsigned char byteRed;

  public:
	inline CImage2DByteRGBPixel() {}
	inline CImage2DByteRGBPixel(unsigned char red, unsigned char green, unsigned char blue)
	{
		byteRed = red;
		byteGreen = green;
		byteBlue = blue;
	}
	inline void Set(unsigned char red, unsigned char green, unsigned char blue)
	{
		byteRed = red;
		byteGreen = green;
		byteBlue = blue;
	}
	inline bool operator ==(const CImage2DByteRGBPixel &rgb) const
	{
		if (byteRed==rgb.byteRed && byteGreen==rgb.byteGreen && byteBlue==rgb.byteBlue)
			return true;
		else
			return false;
	}
	inline bool operator !=(const CImage2DByteRGBPixel &rgb) const
	{
		if (byteRed!=rgb.byteRed || byteGreen!=rgb.byteGreen || byteBlue!=rgb.byteBlue)
			return true;
		else
			return false;
	}
	inline unsigned char ToGrayScale() const {return (unsigned char)(((int)byteRed+(int)byteGreen+(int)byteBlue)/3);}
};

class CImage2DIntRGBPixel
{
  public:
	int iBlue;
	int iGreen;
	int iRed;

  public:
	inline CImage2DIntRGBPixel() {}
	inline CImage2DIntRGBPixel(int red, int green, int blue)
	{
		iRed = red;
		iGreen = green;
		iBlue = blue;
	}
	inline CImage2DIntRGBPixel(const CImage2DByteRGBPixel &byteRGBPixel)
	{
		iRed = (int)byteRGBPixel.byteRed;
		iGreen = (int)byteRGBPixel.byteGreen;
		iBlue = (int)byteRGBPixel.byteBlue;
	}
	inline void Set(int red, int green, int blue)
	{
		iRed = red;
		iGreen = green;
		iBlue = blue;
	}
	inline CImage2DIntRGBPixel operator +(const CImage2DIntRGBPixel &rgbArg) const
	{
		CImage2DIntRGBPixel rgbResult;
		rgbResult.iRed = iRed+rgbArg.iRed;
		rgbResult.iGreen = iGreen+rgbArg.iGreen;
		rgbResult.iBlue = iBlue+rgbArg.iBlue;
		return rgbResult;
	}
	inline CImage2DIntRGBPixel operator -(const CImage2DIntRGBPixel &rgbArg) const
	{
		CImage2DIntRGBPixel rgbResult;
		rgbResult.iRed = iRed-rgbArg.iRed;
		rgbResult.iGreen = iGreen-rgbArg.iGreen;
		rgbResult.iBlue = iBlue-rgbArg.iBlue;
		return rgbResult;
	}
	inline CImage2DIntRGBPixel operator *(float fArg) const
	{
		CImage2DIntRGBPixel rgbResult;
		rgbResult.iRed = (int)(iRed*fArg);
		rgbResult.iGreen = (int)(iGreen*fArg);
		rgbResult.iBlue = (int)(iBlue*fArg);
		return rgbResult;
	}
	inline CImage2DIntRGBPixel operator *(const CImage2DIntRGBPixel &rgbArg) const
	{
		CImage2DIntRGBPixel rgbResult;
		rgbResult.iRed = iRed*rgbArg.iRed;
		rgbResult.iGreen = iGreen*rgbArg.iGreen;
		rgbResult.iBlue = iBlue*rgbArg.iBlue;
		return rgbResult;
	}
	inline CImage2DIntRGBPixel operator /(float fArg) const
	{
		CImage2DIntRGBPixel rgbResult;
		rgbResult.iRed = (int)(iRed/fArg);
		rgbResult.iGreen = (int)(iGreen/fArg);
		rgbResult.iBlue = (int)(iBlue/fArg);
		return rgbResult;
	}
	inline friend CImage2DIntRGBPixel operator *(float fArg, const CImage2DIntRGBPixel &rgbArg) {return rgbArg*fArg;}
	inline friend CImage2DIntRGBPixel operator /(float fArg, const CImage2DIntRGBPixel &rgbArg) {return rgbArg/fArg;}

	inline CImage2DIntRGBPixel &operator +=(const CImage2DIntRGBPixel &rgbArg)
	{
		iRed += rgbArg.iRed;
		iGreen += rgbArg.iGreen;
		iBlue += rgbArg.iBlue;
		return *this;
	}
	inline CImage2DIntRGBPixel &operator -=(const CImage2DIntRGBPixel &rgbArg)
	{
		iRed -= rgbArg.iRed;
		iGreen -= rgbArg.iGreen;
		iBlue -= rgbArg.iBlue;
		return *this;
	}
	inline CImage2DIntRGBPixel &operator *=(float fArg)
	{
		iRed = (int)(iRed*fArg);
		iGreen = (int)(iGreen*fArg);
		iBlue = (int)(iBlue*fArg);
		return *this;
	}
	inline CImage2DIntRGBPixel &operator /=(float fArg)
	{
		iRed = (int)(iRed/fArg);
		iGreen = (int)(iGreen/fArg);
		iBlue = (int)(iBlue/fArg);
		return *this;
	}
	inline CTriplet<int> ToTripletInt() const
	{
		return CTriplet<int>(iRed, iGreen, iBlue);
	}
};

class CImage2DFloatRGBPixel
{
  public:
	float fBlue;
	float fGreen;
	float fRed;

  public:
	inline CImage2DFloatRGBPixel() {}
	inline CImage2DFloatRGBPixel(float red, float green, float blue)
	{
		fRed = red;
		fGreen = green;
		fBlue = blue;
	}
	inline CImage2DFloatRGBPixel(const CImage2DByteRGBPixel &byteRGBPixel)
	{
		fRed = (float)byteRGBPixel.byteRed;
		fGreen = (float)byteRGBPixel.byteGreen;
		fBlue = (float)byteRGBPixel.byteBlue;
	}
	inline CImage2DFloatRGBPixel(const CImage2DIntRGBPixel &iRGBPixel)
	{
		fRed = (float)iRGBPixel.iRed;
		fGreen = (float)iRGBPixel.iGreen;
		fBlue = (float)iRGBPixel.iBlue;
	}

	inline void Set(float red, float green, float blue)
	{
		fRed = red;
		fGreen = green;
		fBlue = blue;
	}

	// Set to zero is used in template class CArray2D (Convolve)
	// Allows to call Convolve() on CArray2D<CImage2DFloatRGBPixel>
	inline CImage2DFloatRGBPixel &operator =(float fValue)
	{
		fRed = fValue;
		fGreen = fValue;
		fBlue = fValue;
		return *this;
	}

	inline float ToGrayScale() const {return (fRed+fGreen+fBlue)/3.0f;}

	inline CImage2DByteRGBPixel ToByteRGBPixel() const
	{
		CImage2DByteRGBPixel byteRGBPixel;
		byteRGBPixel.byteRed   = (unsigned char)std::max(std::min(fRed, 255.0f), 0.0f);
		byteRGBPixel.byteGreen = (unsigned char)std::max(std::min(fGreen, 255.0f), 0.0f);
		byteRGBPixel.byteBlue  = (unsigned char)std::max(std::min(fBlue, 255.0f), 0.0f);
		return byteRGBPixel;
	}
	inline CImage2DIntRGBPixel ToIntRGBPixel() const
	{
		CImage2DIntRGBPixel iRGBPixel;
		iRGBPixel.iRed   = (int)fRed;
		iRGBPixel.iGreen = (int)fGreen;
		iRGBPixel.iBlue  = (int)fBlue;
		return iRGBPixel;
	}
	inline CTriplet<float> ToTripletFloat() const
	{
		return CTriplet<float>(fRed, fGreen, fBlue);
	}
	inline CImage2DFloatRGBPixel operator +(const CImage2DFloatRGBPixel &rgbArg) const
	{
		CImage2DFloatRGBPixel rgbResult;
		rgbResult.fRed = fRed+rgbArg.fRed;
		rgbResult.fGreen = fGreen+rgbArg.fGreen;
		rgbResult.fBlue = fBlue+rgbArg.fBlue;
		return rgbResult;
	}
	inline CImage2DFloatRGBPixel operator -(const CImage2DFloatRGBPixel &rgbArg) const
	{
		CImage2DFloatRGBPixel rgbResult;
		rgbResult.fRed = fRed-rgbArg.fRed;
		rgbResult.fGreen = fGreen-rgbArg.fGreen;
		rgbResult.fBlue = fBlue-rgbArg.fBlue;
		return rgbResult;
	}
	inline CImage2DFloatRGBPixel operator *(float fArg) const
	{
		CImage2DFloatRGBPixel rgbResult;
		rgbResult.fRed = fRed*fArg;
		rgbResult.fGreen = fGreen*fArg;
		rgbResult.fBlue = fBlue*fArg;
		return rgbResult;
	}
	inline CImage2DFloatRGBPixel operator *(const CImage2DFloatRGBPixel &rgbArg) const
	{
		CImage2DFloatRGBPixel rgbResult;
		rgbResult.fRed = fRed*rgbArg.fRed;
		rgbResult.fGreen = fGreen*rgbArg.fGreen;
		rgbResult.fBlue = fBlue*rgbArg.fBlue;
		return rgbResult;
	}
	inline CImage2DFloatRGBPixel operator /(float fArg) const
	{
		CImage2DFloatRGBPixel rgbResult;
		rgbResult.fRed = fRed/fArg;
		rgbResult.fGreen = fGreen/fArg;
		rgbResult.fBlue = fBlue/fArg;
		return rgbResult;
	}
	inline friend CImage2DFloatRGBPixel operator *(float fArg, const CImage2DFloatRGBPixel &rgbArg) {return rgbArg*fArg;}
	inline friend CImage2DFloatRGBPixel operator /(float fArg, const CImage2DFloatRGBPixel &rgbArg) {return rgbArg/fArg;}

	inline CImage2DFloatRGBPixel &operator +=(const CImage2DFloatRGBPixel &rgbArg)
	{
		fRed += rgbArg.fRed;
		fGreen += rgbArg.fGreen;
		fBlue += rgbArg.fBlue;
		return *this;
	}
	inline CImage2DFloatRGBPixel &operator -=(const CImage2DFloatRGBPixel &rgbArg)
	{
		fRed -= rgbArg.fRed;
		fGreen -= rgbArg.fGreen;
		fBlue -= rgbArg.fBlue;
		return *this;
	}
	inline CImage2DFloatRGBPixel &operator *=(float fArg)
	{
		fRed *= fArg;
		fGreen *= fArg;
		fBlue *= fArg;
		return *this;
	}
	inline CImage2DFloatRGBPixel &operator /=(float fArg)
	{
		fRed /= fArg;
		fGreen /= fArg;
		fBlue /= fArg;
		return *this;
	}
	inline friend CImage2DFloatRGBPixel sqrt(const CImage2DFloatRGBPixel &rgbArg)
	{
		CImage2DFloatRGBPixel rgbResult;
		rgbResult.fRed = sqrt(rgbArg.fRed);
		rgbResult.fGreen = sqrt(rgbArg.fGreen);
		rgbResult.fBlue = sqrt(rgbArg.fBlue);
		return rgbResult;
	}
	inline float L1Norm() const
	{
		return fabs(fRed)+fabs(fGreen)+fabs(fBlue);
	}
	inline float L2Norm() const
	{
		return sqrt(fRed*fRed + fGreen*fGreen + fBlue*fBlue);
	}
	inline float InfiniteNorm() const
	{
		return max(fabs(fRed), max(fabs(fGreen), fabs(fBlue)));
	}
};

typedef float (*typePixelDistanceFunction)(const CImage2DFloatRGBPixel &, const CImage2DFloatRGBPixel &);
typedef CImage2DByteRGBPixel (*typeByteRGBPixelTransformFunction)(const CImage2DByteRGBPixel &);



class CImage2D : protected CArray1D<CImage2DPixel>
{
  // Nested type for image format
  public:
	typedef enum {FORMAT_BMP, FORMAT_JPG, FORMAT_PNG} typeFormat;

  // Static members
#ifdef IMAGE2D_RUNTIME_CHECK
  protected:
	// Pixel returned in case of invalid access with function Pixel()
	static CImage2DPixel pixelError;
#endif

  // Members
  protected:
	static CArray1D<unsigned char> vectBitsBitmapInfo;
	int iWidth, iHeight, iBytesPerLine, iBitsPerPixel;

  // Member functions
  public:
	static void StartImage2D();

	inline CImage2D():CArray1D<CImage2DPixel>() {iWidth = iHeight = iBytesPerLine = iBitsPerPixel = 0;}
	CImage2D(const CImage2D &);
	inline ~CImage2D() {Empty();}

	// Get main properties of image (read only)
	inline int GetWidth() const {return iWidth;}
	inline int GetHeight() const {return iHeight;}
	inline int GetBytesPerLine() const {return iBytesPerLine;}
	inline int GetBitsPerPixel() const {return iBitsPerPixel;}
	inline CImage2DPixel *GetBits() const {return pElements;}
	inline CCouple<int> GetSize() const {return CCouple<int>(iWidth, iHeight);}
	inline int GetOffset(int x, int y) const {return y*iBytesPerLine+x*iBitsPerPixel/8;}
	inline int GetOffset(const CCouple<int> &coord) const {return coord.y*iBytesPerLine+coord.x*iBitsPerPixel/8;}

	bool Create(int, int, int);
	bool Load(const char *);
	bool Save(const char *, typeFormat) const;
	void Empty();

	CImage2D &operator =(const CImage2D &);

	// Access single element (read/write)
	inline CImage2DPixel &Pixel(int x, int y) const
	{
		#ifdef IMAGE2D_RUNTIME_CHECK
		if (iBitsPerPixel!=8)
		{
			cerr<<"ERROR in CImage2D::Pixel(int,int): image is not 8-bit"<<endl;
			return pixelError;
		}
		if (x<0 || x>=iWidth || y<0 || y>=iHeight)
		{
			cerr<<"ERROR in CImage2D::Pixel(int, int): accessing pixel ("<<x<<","<<y<<") out of range [0.."<<iWidth-1<<"]x[0.."<<iHeight-1<<"]"<<endl;
			return pixelError;
		}
		#endif
		return pElements[y*iBytesPerLine+x];
	}
	inline CImage2DPixel &Pixel(const CCouple<int> &coord) const
	{
		#ifdef IMAGE2D_RUNTIME_CHECK
		if (iBitsPerPixel!=8)
		{
			cerr<<"ERROR in CImage2D::Pixel(const CCouple<int> &): image is not 8-bit"<<endl;
			return pixelError;
		}
		if (coord.x<0 || coord.x>=iWidth || coord.y<0 || coord.y>=iHeight)
		{
			cerr<<"ERROR in CImage2D::Pixel(const CCouple<int> &): accessing pixel ("<<coord.x<<","<<coord.y<<") out of range [0.."<<iWidth-1<<"]x[0.."<<iHeight-1<<"]"<<endl;
			return pixelError;
		}
		#endif
		return pElements[coord.y*iBytesPerLine+coord.x];
	}

	// Linear interpolation between pixels (read only)
	inline float GetPixelInterpolate(float x, float y) const
	{
		float dx, dy;
		int xi, yi;
		CImage2DPixel *pPixelTemp;
		float pixelInterpolate;

		#ifdef IMAGE2D_RUNTIME_CHECK
		if (iBitsPerPixel!=8)
		{
			cerr<<"ERROR in CImage2D::GetPixelInterpolate(float,float): image is not 8-bit"<<endl;
			return 0.0f;
		}
		if (x<0.0f || x>=(float)(iWidth-1) || y<0.0f || y>=(float)(iHeight-1))
		{
			cerr<<"ERROR in CImage2D::GetPixelInterpolate(float,float): accessing pixel ("<<x<<","<<y<<") out of range [0.."<<iWidth-1<<"[x[0.."<<iHeight-1<<"["<<endl;
			return 0.0f;
		}
		#endif

		xi = (int)x;
		yi = (int)y;
		dx = x-floor(x);
		dy = y-floor(y);

		// Get address of nearest pixel with lower integer coordinates
		pPixelTemp = pElements + yi*iBytesPerLine + xi;

		pixelInterpolate =
			(1.0f-dx)*(1.0f-dy) * (float)pPixelTemp[0] +
			dx*(1.0f-dy)        * (float)pPixelTemp[1] +
			(1.0f-dx)*dy        * (float)pPixelTemp[iBytesPerLine] +
			dx*dy               * (float)pPixelTemp[iBytesPerLine+1];
		return pixelInterpolate;
	}
	inline float GetPixelInterpolate(const CCouple<float> &coord) const {return GetPixelInterpolate(coord.x,coord.y);}

	inline CImage2DByteRGBPixel GetByteRGBPixel(int x, int y) const
	{
		CImage2DByteRGBPixel byteRGB;
		CImage2DPixel *pPixel;

		#ifdef IMAGE2D_RUNTIME_CHECK
		if (iBitsPerPixel!=24 && iBitsPerPixel!=32)
		{
			cerr<<"ERROR in CImage2D::GetByteRGBPixel(int,int): image is not RGB"<<endl;
			return byteRGB;
		}
		if (x<0 || x>=iWidth || y<0 || y>=iHeight)
		{
			cerr<<"ERROR in CImage2D::GetByteRGBPixel(int, int): accessing pixel ("<<x<<","<<y<<") out of range [0.."<<iWidth-1<<"]x[0.."<<iHeight-1<<"]"<<endl;
			return byteRGB;
		}
		#endif

		pPixel = pElements + y*iBytesPerLine + x*iBitsPerPixel/8;
		byteRGB = *((CImage2DByteRGBPixel *)pPixel);

		return byteRGB;
	}
	inline CImage2DByteRGBPixel GetByteRGBPixel(const CCouple<int> &coord) const
	{
		return GetByteRGBPixel(coord.x, coord.y);
	}
	inline void SetByteRGBPixel(int x, int y, const CImage2DByteRGBPixel &byteRGB)
	{
		CImage2DPixel *pPixel;

		#ifdef IMAGE2D_RUNTIME_CHECK
		if (iBitsPerPixel!=24 && iBitsPerPixel!=32)
		{
			cerr<<"ERROR in CImage2D::SetByteRGBPixel(...): image is not RGB"<<endl;
			return;
		}
		if (x<0 || x>=iWidth || y<0 || y>=iHeight)
		{
			cerr<<"ERROR in CImage2D::SetByteRGBPixel(...): accessing pixel ("<<x<<","<<y<<") out of range [0.."<<iWidth-1<<"]x[0.."<<iHeight-1<<"]"<<endl;
			return;
		}
		#endif

		pPixel = pElements + y*iBytesPerLine + x*iBitsPerPixel/8;
		*((CImage2DByteRGBPixel *)pPixel) = byteRGB;
	}
	inline void SetByteRGBPixel(const CCouple<int> &coord, const CImage2DByteRGBPixel &byteRGB)
	{
		SetByteRGBPixel(coord.x, coord.y, byteRGB);
	}

	inline CImage2DIntRGBPixel GetIntRGBPixel(int x, int y) const
	{
		CImage2DIntRGBPixel iRGB;
		CImage2DPixel *pPixel;

		#ifdef IMAGE2D_RUNTIME_CHECK
		if (iBitsPerPixel!=24 && iBitsPerPixel!=32)
		{
			cerr<<"ERROR in CImage2D::GetIntRGBPixel(int,int): image is not RGB"<<endl;
			return iRGB;
		}
		if (x<0 || x>=iWidth || y<0 || y>=iHeight)
		{
			cerr<<"ERROR in CImage2D::GetIntRGBPixel(int,int): accessing pixel ("<<x<<","<<y<<") out of range [0.."<<iWidth-1<<"]x[0.."<<iHeight-1<<"]"<<endl;
			return iRGB;
		}
		#endif

		pPixel = pElements + y*iBytesPerLine + x*iBitsPerPixel/8;
		iRGB.iBlue = (int)pPixel[0];
		iRGB.iGreen = (int)pPixel[1];
		iRGB.iRed = (int)pPixel[2];

		return iRGB;
	}
	inline CImage2DIntRGBPixel GetIntRGBPixel(const CCouple<int> &coord) const
	{
		return GetIntRGBPixel(coord.x, coord.y);
	}
	inline CImage2DFloatRGBPixel GetFloatRGBPixel(int x, int y) const
	{
		CImage2DFloatRGBPixel fRGB;
		CImage2DPixel *pPixel;

		#ifdef IMAGE2D_RUNTIME_CHECK
		if (iBitsPerPixel!=24 && iBitsPerPixel!=32)
		{
			cerr<<"ERROR in CImage2D::GetFloatRGBPixel(int,int): image is not RGB"<<endl;
			return fRGB;
		}
		if (x<0 || x>=iWidth || y<0 || y>=iHeight)
		{
			cerr<<"ERROR in CImage2D::GetFloatRGBPixel(int,int): accessing pixel ("<<x<<","<<y<<") out of range [0.."<<iWidth-1<<"]x[0.."<<iHeight-1<<"]"<<endl;
			return fRGB;
		}
		#endif

		pPixel = pElements + y*iBytesPerLine+x*iBitsPerPixel/8;
		fRGB.fBlue = (float)pPixel[0];
		fRGB.fGreen = (float)pPixel[1];
		fRGB.fRed = (float)pPixel[2];

		return fRGB;
	}
	inline CImage2DFloatRGBPixel GetFloatRGBPixel(const CCouple<int> &coord) const
	{
		return GetFloatRGBPixel(coord.x, coord.y);
	}

	// Linear interpolation between RGB pixels
	inline CImage2DFloatRGBPixel GetRGBPixelInterpolate(float x, float y) const
	{
		CImage2DFloatRGBPixel fRGB;
		float dx, dy;
		int xi, yi;

		#ifdef IMAGE2D_RUNTIME_CHECK
		if (iBitsPerPixel!=24 && iBitsPerPixel!=32)
		{
			cerr<<"ERROR in CImage2D::GetRGBPixelInterpolate(float,float): image is not RGB"<<endl;
			return fRGB;
		}
		if (x<0.0f || x>=(float)(iWidth-1) || y<0.0f || y>=(float)(iHeight-1))
		{
			cerr<<"ERROR in CImage2D::GetRGBPixelInterpolate(float,float): accessing pixel ("<<x<<","<<y<<") out of range [0.."<<iWidth-1<<"[x[0.."<<iHeight-1<<"["<<endl;
			return fRGB;
		}
		#endif

		xi = (int)x;
		yi = (int)y;
		dx = x-floor(x);
		dy = y-floor(y);

		fRGB =
			(1.0f-dx)*(1.0f-dy) * GetFloatRGBPixel(xi,yi) +
			dx*(1.0f-dy)        * GetFloatRGBPixel(xi+1,yi) +
			(1.0f-dx)*dy        * GetFloatRGBPixel(xi,yi+1) +
			dx*dy               * GetFloatRGBPixel(xi+1,yi+1);
		return fRGB;
	}
	inline CImage2DFloatRGBPixel GetRGBPixelInterpolate(const CCouple<float> &coord) const
	{
		return GetRGBPixelInterpolate(coord.x,coord.y);
	}

	// Draw line in 8-bit image
	// Params: first endpoint, second endpoint, intensity
	void DrawLine(const CCouple<int> &, const CCouple<int> &, CImage2DPixel);

	// Draw filled rectangle in 8-bit image
	// Params: top-left corner, bottom-right corner, intensity
	void DrawFilledRectangle(const CCouple<int> &, const CCouple<int> &, CImage2DPixel col);

	// Draw filled circle in 8-bit image
	// Params: center, radius, intensity
	void DrawFilledCircle(const CCouple<int> &, int, CImage2DPixel col);

	// Draw line in RGB image
	// Params: first endpoint, second endpoint, color
	void DrawLineRGB(const CCouple<int> &, const CCouple<int> &, const CImage2DByteRGBPixel &);

	// Draw filled rectangle in RGB image
	// Params: top-left corner, bottom-right corner, color
	void DrawFilledRectangleRGB(const CCouple<int> &, const CCouple<int> &, const CImage2DByteRGBPixel &);

	// Draw filled circle in RGB image
	// Params: center, radius, color
	void DrawFilledCircleRGB(const CCouple<int> &, int, const CImage2DByteRGBPixel &);

	// Return an 8-bit resized image
	// No interpolation between pixels -> resampling is performed with truncation
	// Params: new width, new height
	CImage2D Resize(int, int) const;

	// Return a RGB resized image
	// No interpolation between pixels -> resampling is performed with truncation
	// Params: new width, new height
	CImage2D ResizeRGB(int, int) const;

	// Return an 8-bit image from a RGB image
	// Intensity is just taken as the average of RGB components
	CImage2D GrayScale() const;

	// Test if an RGB image, stored in 24 or 32 bits, contains only grayscale pixels
	bool IsGrayScaleRGB() const;

	// Convert 8-bit to RGB
	// Param: output RGB image
	void Convert8bitsto24bits(CImage2D &dest24bits) const;

	void Clear(CImage2DPixel);
	void ClearRGB(const CImage2DByteRGBPixel &);

	void Flip();

	void TransformWithLUT(const CArray1D<unsigned char> &);
	void TransformWithFunctionRGB(typeByteRGBPixelTransformFunction);

    // Pixelwise operations on 8-bit images
    CImage2D AbsDiff(const CImage2D &) const;
    CImage2D Min(const CImage2D &) const;
    CImage2D Max(const CImage2D &) const;

    // Flood fill from seed pixel with given intensity
	// Params: initial seed pixel, intensity
	void FloodFill(const CCouple<int> &, CImage2DPixel);

	// Create image from array of real values
	bool CreateFromArray2DFloat(const CArray2D<float> &);
	bool CreateFromArray2DTripletFloatRGB(CArray2D<CTriplet<float> > &);

	// Convert image to array of real values
	bool ConvertToArray2DFloat(CArray2D<float> &) const;
	bool ConvertToArray2DTripletFloatRGB(CArray2D<CTriplet<float> > &) const;
	bool ConvertToArray2DTripletFloatYUV(CArray2D<CTriplet<float> > &) const;
	bool ConvertToArray2DTripletFloatLAB(CArray2D<CTriplet<float> > &) const;
	bool ConvertToArray2DFloatRGBPixel(CArray2D<CImage2DFloatRGBPixel> &) const;

  // Decoders and encoders
  protected:

	// BMP encoder/decoder
	bool DecodeBMP(FILE *);
	bool EncodeBMP(FILE *) const;

	// JPEG encoder/decoder
	#ifdef IMAGE2D_SUPPORT_JPG
	bool DecodeJPG(FILE *);
	bool EncodeJPG(FILE *) const;
	#endif // IMAGE2D_SUPPORT_JPG

	#ifdef IMAGE2D_SUPPORT_PNG
	bool DecodePNG(FILE *);
	bool EncodePNG(FILE *) const;
	#endif // IMAGE2D_SUPPORT_PNG
};

#endif
