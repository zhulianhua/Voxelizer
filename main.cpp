#include <iostream>
#include <cstdlib>
#include "Voxelizer.h"

int main(int argc, char* argv[])
 {
     // usage
     if (argc != 2) {
         cerr << "Usage voxelizer.exe <filename>.stl" << endl;
         exit(EXIT_FAILURE);
     }
    string fname(argv[1]);

    // import geometry
    Geometry geo(fname);
    // The geometry bound, min/max corners coordinates will be reported.
    // These information will be used to determine the appropriate bound the grid
    // For this example we use the cube.stl file
    // Tese reported information are
    //
    //GEO bound = Bbox(
    //minCorner:    (-1.0, 1.0, 1.0)
    //maxCorner:    (-1.0, 1.0, 1.0)
    //Extend:       ( 2,   2,   2  )
    //)
    // Now we need create the bound of the **grid**. 
    // We set it larger then the bound of the **geometry**.
    // Let's say: 
    // minCorner (-1,-1,-0.5)
    // maxCorner ( 1, 1, 1.5)
    V3 gridMinCorner(-2.0, -2.0, -2.0);
    V3 gridMaxCorner( 2.0,  2.0,  2.0);
    V3 gridExtend = gridMaxCorner - gridMinCorner;

    // and we want a 33 voxels along X direction,
    // so dx = (1-(-1))/32;
    double dx = gridExtend.x / 33;
    // and the grid number in each direction can be determined
    // Why not 32??? The voxillization implemntation is not robust enough.
    // If use 32, then the polar points causes wrong vexillization in the sphere axis

    // set grid bounding corner
    //int3 gridNum = int3(gridExtend.x / dx, gridExtend.y / dx, gridExtend.z /dx);

    // create the grid finally
    GridBox grid(gridMinCorner, dx, gridNum);

    // generate voxilzer from geometry and gridbox
    Voxelizer vox(geo, grid);

    // get flag
    const char* flag = vox.get_flag();
    gridNum = grid.get_gridNum();

    // count solid flag
    int count = 0;
    for (int i = 0; i < gridNum.nx * gridNum.ny * gridNum.nz; i++)
    {
        //if (*(flag+i) == 0) count++;
        if (flag[i] == 1) count++;
    }
    cout << "count of flag == 1 : " << count << endl;

    // wirte vtk file of flag, use paraview to view the flag data
    vox.write_vtk_image();

    return 0;
}
