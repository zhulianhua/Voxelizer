// voxelizer.cpp : Defines the entry point for the console application.
//

#include <string>
#include <fstream>
#include <iostream>
#include <memory>
#include <vector>
#include "stdafx.h"
#include "Voxelizer.h"

using namespace std;

Voxelizer::Voxelizer(Geometry& geo_, GridBox& grid_) : geo(geo_), grid(grid_) {
	int3 gridNum = grid.get_gridNum();
	flag = new char[gridNum.nx * gridNum.ny * gridNum.nz];

	for (int iz = 0; iz < gridNum.nz; iz++)
	{
		vector<Triangle>& tris = get_relevant_triangles(iz);
		vector<Line>& lines = get_z_sections(iz, tris);
		for (int iy = 0; iy < gridNum.ny; iy++)
		{
			vector<int>& xids = get_xid_cross(iy, lines);
		}
	}
}

void Voxelizer::read_stl_file(string fname) {
	ifstream stlFile(fname.c_str(), ios::in | ios::binary);
	char headInfo[80] = "";
	// read 80 byte header
	if (stlFile)
	{
		stlFile.read(headInfo, 80);
		cout << ">> Reading " << headInfo << endl;
		char nTriRaw[4];
		stlFile.read(nTriRaw, 4);
		unsigned numTri = *((unsigned*)nTriRaw);
		triangles.resize(numTri);
		for (auto&& tri : triangles)
		{
			char triRaw[50];
			stlFile.read(triRaw, 50);
			V3 norm(triRaw);
			V3 p1(triRaw+12);
			V3 p2(triRaw+24);
			V3 p3(triRaw+36);
			tri = Triangle(norm, p1, p2, p3);
		}
	}
	else
	{
		cerr << "File openning error!" << endl;
	}
}

void Voxelizer::scale_shift_mesh(double scale_, V3 shift_) {
	// keep a track of the current shift and scale v.s. to the STL file
	shift += shift_;
	scale *= scale_;
	for (auto&& tri : triangles)
	{
		tri.p1 += shift_;
		tri.p1 *= scale_;
		tri.p2 += shift_;
		tri.p2 *= scale_;
		tri.p3 += shift_;
		tri.p3 *= scale_;
		tri.norm *= scale_;
	}
	minCorner += shift_;
	minCorner *= scale_;
	maxCorner += shift_;
	maxCorner *= scale_;
}

void Voxelizer::set_min_max_corner() {
	double minX, minY, minZ, maxX, maxY, maxZ;
	minX = minY = minZ =  1.0e30;
	maxX = maxY = maxZ = -1.0e30;
	for (auto&& tri : triangles)
	{
		minX = (tri.p1.x < minX) ? tri.p1.x : minX;
		minX = (tri.p2.x < minX) ? tri.p2.x : minX;
		minX = (tri.p3.x < minX) ? tri.p3.x : minX;
		minY = (tri.p1.y < minY) ? tri.p1.y : minY;
		minY = (tri.p2.y < minY) ? tri.p2.y : minY;
		minY = (tri.p3.y < minY) ? tri.p3.y : minY;
		minZ = (tri.p1.z < minZ) ? tri.p1.z : minZ;
		minZ = (tri.p2.z < minZ) ? tri.p2.z : minZ;
		minZ = (tri.p3.z < minZ) ? tri.p3.z : minZ;

		maxX = (tri.p1.x > maxX) ? tri.p1.x : maxX;
		maxX = (tri.p2.x > maxX) ? tri.p2.x : maxX;
		maxX = (tri.p3.x > maxX) ? tri.p3.x : maxX;
		maxY = (tri.p1.y > maxY) ? tri.p1.y : maxY;
		maxY = (tri.p2.y > maxY) ? tri.p2.y : maxY;
		maxY = (tri.p3.y > maxY) ? tri.p3.y : maxY;
		maxZ = (tri.p1.z > maxZ) ? tri.p1.z : maxZ;
		maxZ = (tri.p2.z > maxZ) ? tri.p2.z : maxZ;
		maxZ = (tri.p3.z > maxZ) ? tri.p3.z : maxZ;
	}
	minCorner = V3(minX, minY, minZ);
	maxCorner = V3(maxX, maxY, maxZ);
}

vector<Triangle>& Voxelizer::get_relevant_triangles(int iz_) const
{
}


int main(int argc, char* argv[])
 {
	string fname(argv[1]);
	Voxelizer vox(fname, 100);
	cout << "Num of triangle " << vox.triangles.size() << endl;
	Triangle tri = vox.triangles[0];
	//for (auto&& tri : vox.triangles)
	//	cout << tri << endl;
	cout.precision(10);
	cout << "0th tri : " << tri << endl;
	cout << "Max corner " <<  vox.maxCorner << endl;
	cout << "Min corner " <<  vox.minCorner << endl;
	cout << "Extend  = " << V3(vox.maxCorner - vox.minCorner) << endl;
    return 0;
}