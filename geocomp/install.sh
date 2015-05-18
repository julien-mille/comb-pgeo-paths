#!/bin/bash

if [ ! -d "build" ]; then
    mkdir build;
fi
cd build;
echo "Run cmake";
cmake ..;
echo "Run make";
make;
