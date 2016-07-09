// voxelizer.cpp : Defines the entry point for the console application.
//

#include <string>
#include <fstream>
#include <iostream>
#include <memory>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include "stdafx.h"
#include "Voxelizer.h"

using namespace std;

// global functions
// determin z plane / triangle relative position type
vector<int> tri_plane(Triangle& tri_, double z_){
	vector<int> ans;
	if (tri_.p1.z > z_) ans.push_back(1);
	if (tri_.p2.z > z_) ans.push_back(2);
	if (tri_.p3.z > z_) ans.push_back(3);
	return ans;
}

// get the section line of triangle with a z plane
Line tri_sec_plane(Triangle& tri_, vector<int>& aboveP_, double z_) {
	// find which points are above
	//if (tri_plane(tri_, z_) == 1);
	//if (ans.size() )
	Line line1, line2;
	if (aboveP_.size() == 1) {
		if (aboveP_[0] == 1) {
			line1 = Line(tri_.p1, tri_.p2);
			line2 = Line(tri_.p1, tri_.p3);
		}
		else if (aboveP_[0] == 2) {
			line1 = Line(tri_.p2, tri_.p1);
			line2 = Line(tri_.p3, tri_.p3);
		}
		else {
			line1 = Line(tri_.p3, tri_.p1);
			line2 = Line(tri_.p3, tri_.p2);
		}
	}
	else {  // 2 above points
		// find the below point first
		if (aboveP_[0] == 1 && aboveP_[1] == 2) {
			line1 = Line(tri_.p3, tri_.p1);
			line2 = Line(tri_.p3, tri_.p2);
		} 
		else if (aboveP_[0] == 1 && aboveP_[1] == 3) {
			line1 = Line(tri_.p2, tri_.p3);
			line2 = Line(tri_.p2, tri_.p1);
		}
		else {
			line1 = Line(tri_.p1, tri_.p2);
			line2 = Line(tri_.p1, tri_.p3);
		}
	}
	//
	V3 p1, p2;
	p1 = line1.p_cross_z_plane(z_);
	p2 = line2.p_cross_z_plane(z_);
	return Line(p1, p2);
}


ostream& operator<<(ostream& os_, Bbox& box_)
{
	os_ << "Bbox(" << endl 
		<< "minCorner: \t" << box_.minCorner << endl
		<< "maxCorner: \t" << box_.maxCorner << endl
		<< "Extend: \t" << box_.get_extend() << endl << ")" << endl;
		
	return os_;
}

Line::Line(V3& p1_, V3& p2_) : p1(p1_), p2(p2_) {};

V3 Line::p_cross_z_plane(double z) const {
	V3 ans = p1 + (p2 - p1) * (z - p1.z) / (p2.z - p1.z);
	return ans;
}

V3 Line::p_cross_y_plane(double y) const {
	V3 ans = p1 + (p2 - p1) * (y - p1.y) / (p2.y - p1.y);
	return ans;
}

V3 Line::p_cross_x_plane(double x) const {
	V3 ans = p1 + (p2 - p1) * (x - p1.x) / (p2.x - p1.x);
	return ans;
}

Geometry::Geometry(string fname) {
	read_stl_file(fname);
	cout << "Num of triangle " << get_num_tri() << endl;
	set_bound();
}

Geometry::~Geometry() {}

GridBox::GridBox(V3& minCorner_, V3& maxCorner_, double dx_) 
	: Bbox(minCorner_, maxCorner_), dx(dx_)
{
	// the size in x/y/z direction has to be multiple of dx
	V3 extend = maxCorner - minCorner;
	int nxI = (int)(extend.x / dx);
	int nyI = (int)(extend.y / dx);
	int nzI = (int)(extend.z / dx);
	double xRes = extend.x - nxI * dx;
	double yRes = extend.y - nyI * dx;
	double zRes = extend.z - nzI * dx;
	if (xRes > 1e-6 * dx || yRes > 1e-6 * dx || zRes > 1e-6 * dx)
		exit(EXIT_FAILURE);
	gridNum = int3(nxI, nyI, nzI);
}

GridBox::GridBox(V3& minCorner_, double dx_, int3 gridNum_)
	: dx(dx_), gridNum(gridNum_)
{
	double lx = dx_ * gridNum.nx;
	double ly = dx_ * gridNum.ny;
	double lz = dx_ * gridNum.nz;
	minCorner = minCorner_;
	maxCorner = minCorner + V3(lx, ly, lz);
}

void Geometry::read_stl_file(string fname) {
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

void Geometry::scale_shift(double scale_, V3 shift_) {
	// keep a track of the current shift and scale v.s. to the STL file
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
	// min/max corner of the bound
	V3 minCorner = bound.get_minCorner();
	V3 maxCorner = bound.get_maxCorner();
	minCorner += shift_;
	minCorner *= scale_;
	bound = Bbox(minCorner, maxCorner);
}

void Geometry::set_bound() {
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
	V3 minCorner = V3(minX, minY, minZ);
	V3 maxCorner = V3(maxX, maxY, maxZ);
	bound = Bbox(minCorner, maxCorner);
}

Voxelizer::Voxelizer(Geometry& geo_, GridBox& grid_) : geo(geo_), grid(grid_) {
	// scale (in the z direction) and shift the geometry to fit the grid
	Bbox bound = geo.get_bound();
	double scale = grid.get_extend().z / bound.get_extend().z;
	geo.scale_shift(scale, V3(0.0, 0.0, 0.0));
	V3 shift = grid.get_minCorner() - bound.get_minCorner();
	geo.scale_shift(1, shift);

	int3 gridNum = grid.get_gridNum();
	int nx = gridNum.nx;
	int ny = gridNum.ny;
	int nz = gridNum.nz;
	long numTotal = nx * ny * nz;
	flag = new char[numTotal];
	for (int i = 0; i < numTotal; i++)
		flag[i] = 0;

	for (int iz = 1; iz <= nz; iz++)
	{
		vector<Line> lines;
		get_z_sections(lines, iz);
		for (int iy = 1; iy <= ny; iy++)
		{
			vector<int> xids;
			get_xid_cross(xids, iy, lines);
			int n_change = 0;
			bool isBlack = true;
			for (int ix = 0; ix < nx; ix++)
				if (ix == xids[n_change]) {
					if (isBlack) {
						isBlack = !isBlack;
						flag[(iz - 1) * nx * ny + (iy - 1) * nx + ix] = isBlack;
					}
					else {
						flag[(iz - 1) * nx * ny + (iy - 1) * nx + ix] = isBlack;
						isBlack = !isBlack;
					}
					n_change++;
				}
				else{
					flag[(iz - 1) * nx * ny + (iy - 1) * nx + ix] = isBlack;
				}
		}
	}
}

Voxelizer::~Voxelizer()
{
	cout << "Deleting the flag memory!" << endl;
	delete []flag;
}

void Voxelizer::get_relevant_triangles(vector<Triangle>& tri_, int iz_) const {
	// case 0: 3 points below, 0 points above
	// case 1: 2 points below, 1 points above
	// case 2: 1 points below, 2 points above
	// case 3: 0 points below, 3 points above
	for (int i = 0; i < geo.get_num_tri(); i++) {
		Triangle triTmp = geo.get_tri(i);
		double dx = grid.get_dx();
		double zmin = grid.get_minCorner().z;
		vector<int> aboveP = tri_plane(triTmp, iz_ * dx + zmin);
		if (aboveP.size() > 0 && aboveP.size() < 3) {
			tri_.push_back(triTmp);
		}
	}
}

void Voxelizer::get_z_sections(vector<Line>& lines_, int iz_, vector<Triangle>& tris_) const {}

void Voxelizer::get_z_sections(vector<Line>& lines_, int iz_) const {
	lines_.clear();
	for (int i = 0; i < geo.get_num_tri(); i++) {
		Triangle triTmp = geo.get_tri(i);
		double z = grid.get_dx() * iz_ + grid.get_minCorner().z;
		vector<int> aboveP = tri_plane(triTmp, z);
		if (aboveP.size() > 0 && aboveP.size() < 3) {
			Line lineTmp = tri_sec_plane(triTmp, aboveP, z);
			lines_.push_back(lineTmp);
		}
	}
}

void Voxelizer::get_xid_cross(vector<int>& xids_, int iy_, vector<Line>& lines_) const {
	V3 minCorner = grid.get_minCorner();
	double dx = grid.get_dx();
	double y = iy_ * dx + minCorner.y;
	xids_.clear();
	for (auto && lineTmp : lines_) 
		if ((lineTmp.p1.y - y) * (lineTmp.p2.y - y) < 0.0) { // line cross y
			V3 cross = lineTmp.p_cross_y_plane(y);
			xids_.push_back((int)(cross.x - minCorner.x) / dx);
		}
	// sort
	std::sort(xids_.begin(), xids_.end());
}

int main(int argc, char* argv[])
 {
	 cout << argc << endl;
	 if (argc != 2) {
		 cerr << "Usage voxelizer.exe <filename>.stl" << endl;
		 exit(EXIT_FAILURE);
	 }
	string fname(argv[1]);
	// import geometry
	Geometry geo(fname);
	cout << geo.get_bound() << endl;

	// set bounding box
	V3 minCorner(V3::zero);
	V3 maxCorner(1.0, 1.0, 1.0);
	GridBox grid(minCorner, maxCorner, 0.01); // dx = 0.01, so 100 cells in each direction
	cout << "dx = " << grid.get_dx() << endl;

	Voxelizer vox(geo, grid);
    return 0;
}