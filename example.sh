#!/bin/sh

# An example to build & run an example

rm -rf build
mkdir build && cd build
cmake ..
make
./voxelizer ../sphere.stl
