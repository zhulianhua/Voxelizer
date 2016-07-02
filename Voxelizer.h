#ifndef VOXELIZER_H
#define VOXELIZER_H

#include <vector>
#include <string>
#include <memory>
#include "V3.h"
#include "Triangle.h"

typedef struct int3 { int nx,ny, nz;} int3;


class Line {
public:
	Line(V3& p1_, V3& p2_);
	V3 get_length() const;
	V3 p_cross_x_plane(double x) const;
	V3 p_cross_y_plane(double y) const;
	V3 p_cross_z_plane(double z) const;
	V3 p1;
	V3 p2;
};

bool corss_y(Line line_, double y)
{
	return (line_.p1.y - y) * (line_.p2.y - y) > 0.0;
}

class Bbox {
public:
	Bbox(V3& minCorner_, V3& maxCorner_);
	V3 get_minCorner() const;
	V3 get_maxCorner() const;
	V3 get_extend() const;
	~Bbox();
private:
	V3 minCorner;
	V3 maxCorner;
};

class Geometry {
public:
	Geometry(string fname_);
	int get_num_tri() const;
	V3 get_refPoint() const;
	V3 get_bound() const;
	void set_refPoint(V3& refPoint_);
	void scale_shift(double scale_, V3 shift_);
	~Geometry();
private:
	void read_stl_file(string fname);
	V3 refPoint;
	vector<Triangle> triangles;
	Bbox bound;
};

class GridBox : public Bbox {
public:
	GridBox(V3& minCorner_, V3& macCorner_, double dx_);
	GridBox(Bbox& bound_, double dx_);
	GridBox(V3& minCorner_, double dx_, int3 gridNum_);
	double get_dx() const;
	int3 get_gridNum() const;
private:
	double dx;
	int3 gridNum;
};

class Voxelizer {
public:	
	// constructor do the vexelization
	Voxelizer(Geometry& geo_, GridBox& grid_);
	// return the voxelized flag;
	shared_ptr<char> get_flag();
	// ~descruction
private:
	// get triangles that intersect with plane iz_
	vector<Triangle>& get_relevant_triangles(int iz_) const;

	// get intersection of the geometry with the plane iz_;
	vector<Line>& get_z_sections(int iz_, vector<Triangle>& tris_) const;

	// get lines that intersect with line iy_
	vector<Line>& get_revelant_lines(int iy_ ) const;

	// get xid of the cells that intersection with lines at iy_
	vector<int>& get_xid_cross(int iy_, vector<Line>& lines_) const;
	// flag data
	char* flag;

	// geometry
	Geometry geo;
	// grid box
	GridBox grid;
};
#endif VOXELIZER_H
