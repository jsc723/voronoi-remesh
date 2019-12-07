#ifndef VERTEX_H
#define VERTEX_H
#include "config.h"
#include "vector3.h"
#include <vector>
#include <set>

class Vertex
{
public:
	Vertex(void);
	Vertex(Vec3d);
	Vertex(double,double,double);
	~Vertex(void);

	int id;
	Vec3d pos;
	std::set<int> neighborVertexIds;
	void addNeighborVertex(int);
	void delNeighborVertex(int);
	bool hasNeighborVertex(int);
};
#endif 
