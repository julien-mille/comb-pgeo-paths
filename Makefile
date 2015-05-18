#
# Makefile for CombinationOfPiecewiseGeodesicPaths
# Compiler: gcc/MinGW
#
CC = g++
CPPFLAGS = -DNDEBUG -O3 -msse2
LD = g++
LDFLAGS = -DNDEBUG -O3 -msse2
SRCS = arrayndfloat.cpp histogram.cpp image2d.cpp main.cpp path.cpp piecewisegeodesiccombination.cpp voronoigraph.cpp
FINAL_TARGET = combpaths
INCLUDE_DIR = -I./external/include
LIB_DIR = #-L./external/lib/gcc # set directory to ./external/lib/mingw if using MinGW
LIBS = -lpng -lz -ljpeg 

default: $(FINAL_TARGET)

$(FINAL_TARGET): $(SRCS:%.cpp=%.o)
	$(LD) $+ -o $@ $(LDFLAGS) $(LIB_DIR) $(LIBS)

%.o: %.cpp
	$(CC) -c $(CPPFLAGS) $(INCLUDE_DIR) $< -o $@

clean:
	rm -f *.o $(FINAL_TARGET)

