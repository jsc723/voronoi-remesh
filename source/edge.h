#ifndef EDGE_H
#define EDGE_H
#include "vector3.h"
#include <cstdlib>
#include <cstdio>
// edge class
class Edge
{
public:
	Edge(int _u = -99,int _v = -99);
	~Edge(void);
	int id; 
	int v1,v2; // the id of the two vertices
	Vec3d v; // point after collasping
	double costV; // cost of collasping 
};

#endif