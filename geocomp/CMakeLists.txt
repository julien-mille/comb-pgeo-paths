cmake_minimum_required(VERSION 2.6)

project(GeoComp)

set(CMAKE_VERBOSE_MAKEFILE OFF)

set(CMAKE_C_FLAGS_DEBUG "")
set(CMAKE_CXX_FLAGS_DEBUG "")

set(CMAKE_C_FLAGS_RELEASE "")
set(CMAKE_CXX_FLAGS_RELEASE "")

# Set to Release by default
set(CMAKE_BUILD_TYPE Release CACHE INTERNAL "Build type determining compiler flags (Release or Debug)" FORCE)

# The differents option of CImg (independent from platform)
list(APPEND GEOCOMP_CIMG_INCPATH "/home/sbougleux/Documents/code/CImg/")

if (NOT WIN32)
  #set(X11_LIB_PATH "/usr/lib" CACHE PATH "Path to X11")
  option(GEOCOMP_CIMG_X11 "CImg: Use X11 libraries" ON)
  option(GEOCOMP_CIMG_XSHM "CImg: Use XSHM for fast drawing" OFF)
  option(GEOCOMP_CIMG_XRANDR "CImg: Use XRandR for fast drawing" OFF)
  option(GEOCOMP_CIMG_VT100 "CImg: Use vt100 proto for console printing" OFF)
  option(GEOCOMP_CIMG_LAPACK "CImg: Use Lapack library" OFF)
  option(GEOCOMP_MACPORT "Use Macport installation" OFF)
  list(APPEND GEOCOMP_CIMG_LIBS "pthread")
else (NOT WIN32)
  add_definitions("-mwindows")
  list(APPEND GEOCOMP_CIMG_LIBS "gdi32")
endif (NOT WIN32)

# Manage the option and construct CFLAGS and LDFLAGS
if (GEOCOMP_CIMG_X11)
  list(APPEND GEOCOMP_CIMG_LIBS "X11")
  #list(APPEND GEOCOMP_CIMG_LIBSPATH "${X11_LIB_PATH}/lib")
  #list(APPEND GEOCOMP_CIMG_INCPATH "${X11_LIB_PATH}/include")
endif (GEOCOMP_CIMG_X11)

if (GEOCOMP_CIMG_XSHM)
  list(APPEND GEOCOMP_CIMG_FLAGS "-Dcimg_use_xshm")
  list(APPEND GEOCOMP_CIMG_LIBS "Xext")
endif (GEOCOMP_CIMG_XSHM)

if (GEOCOMP_CIMG_XRANDR)
  list(APPEND GEOCOMP_CIMG_FLAGS "-Dcimg_use_xrandr")
  list(APPEND GEOCOMP_CIMG_LIBS "Xrandr")
endif (GEOCOMP_CIMG_XRANDR)

if (GEOCOMP_CIMG_VT100)
  list(APPEND GEOCOMP_CIMG_FLAGS "-Dcimg_use_vt100")
endif (GEOCOMP_CIMG_VT100)

if (GEOCOMP_CIMG_LAPACK)
  list(APPEND GEOCOMP_CIMG_FLAGS "-Dcimg_use_lapack")
  list(APPEND GEOCOMP_CIMG_LIBS "blas" "lapack")  
endif(GEOCOMP_CIMG_LAPACK)

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
  add_definitions("-O0 -g -ansi -Wall -Wextra -pedantic")
else (CMAKE_BUILD_TYPE STREQUAL "Debug")
  #add_definitions("-O3 -ansi -Wall -Wextra -pedantic")
endif (CMAKE_BUILD_TYPE STREQUAL "Debug")

# We now need these library so add it to CIMG_LIBS
#list(APPEND GEOCOMP_CIMG_LIBS "gsl" "gslcblas")

# This is specific for some Mac configuration
if (GEOCOMP_MACPORT)
  list(APPEND GEOCOMP_CIMG_LIBSPATH "/opt/local/lib")
endif(GEOCOMP_MACPORT)

# The path to the source code
add_subdirectory(example-cimg)
