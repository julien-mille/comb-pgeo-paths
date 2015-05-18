# comb-pgeo-paths
Combination of piecewise-geodesic paths for interactive segmentation

=========================================================================
README for C++ library CombinationOfPiecewiseGeodesicPaths
=========================================================================
Copyright 2015 Julien Mille
=========================================================================
This file is part of CombinationOfPiecewiseGeodesicPaths.

CombinationOfPiecewiseGeodesicPaths is free software: you can redistribute 
it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation, either version 3 of
the License, or (at your option) any later version.

CombinationOfPiecewiseGeodesicPaths is distributed in the hope that it
will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU General Public License
along with GeoComp. If not, see <http://www.gnu.org/licenses/>.
=========================================================================

=========================================================================
HOW TO BUILD AND RUN
=========================================================================

* LINUX

CombinationOfPiecewiseGeodesicPaths was tested under Ubuntu 14.04

If you have libpng, libjpeg and zlib alreadey installed in the default
library directory /usr/lib, then edit image2d.h and to the following
replacements:

#include <png/png.h> -> #include <png.h>
#include <jpeg/jepglib.h> -> #include <jpeglib.h>

Edit the Makefile and comment lines defining INCLUDE_DIR and LIB_DIR

> make
> ./combpaths

Otherwise, if you're using a 32-bit Linux, then leave image2d.h and
Makefile as provided. Run the previous commands. If you're using a 64-bit
and do not wish to install the image libraries, then you should download
and build them on your own and set the Makefile accordingly


* WINDOWS mingw

CombinationOfPiecewiseGeodesicPaths was tested under Windows 7 64bits
with the 32-bit MinGW version shipped with Code::Blocks 13.02

If you're using the 32-bit version of MinGW, then just run the Makefile
and it should be fine

> mingw32-make /f Makefile.mingw
> ./compaths.exe

If you're using a 64-bit version of MinGW, then you'll have to build
libjpeg, libpng and zlib on your own and set Makefile.mingw accordingly


* WINDOWS Microsoft Visual C++

CombinationOfPiecewiseGeodesicPaths was tested under Windows 7 64bits
with Microsoft Visual C++ 2013 Express 32-bit

You may either create a project from scratch and add the existing source
files (select a Console Application if this is the case), or use
Makefile.msvc. Run the MSVC console and enter:

> nmake /F Makefile.msvc
> ./combpaths.exe
