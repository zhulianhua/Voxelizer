#include <iostream>
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

	V3 geoExt = geo.get_bound().get_extend();
	double dx = geoExt.z / 500;
	// set grid bounding corner
	V3 minCorner(0.0, 0.0, 0.0);

	int3 gridNum = int3(geoExt.x / dx + 4, geoExt.y / dx + 4, geoExt.z /dx + 4);
	GridBox grid(minCorner, dx, gridNum);

	//GridBox grid(minCorner, maxCorner, geoExt.z / 100); // dx = 0.01, so 100 cells in each direction

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
