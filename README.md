### STL file to voxel converter written in C++

* See main.cpp for usage
* use cmake to compile
* read only binary stl file

#### TODO
* Improve precision : dectecting on-plane points/lines/triangles


### Example 
The script to run the example is `example.sh`. the grid setting in main.cpp is for this example.
The description of the settings are in the comments in the main.cpp.

`sphere.stl` : A very simple binary stl file, just a sphere centers at (0,0,0) with a radius 1.0:
```
minCorner:  (-1.0, -1.0, -1.0)
maxCorner:  ( 1.0,  1.0,  1.0)
Extend:     ( 2.0,  2.0,  2.0)
```
This information can be reported when first run `./voxelizer ../sphere.stl` in the build directory.
