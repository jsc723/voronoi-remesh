#include "vertex.h"


Vertex::Vertex(void)
{
	id = -99;
	neighborVertexIds.clear();
}

Vertex::Vertex(Vec3d t)
{
	id = -98;
	neighborVertexIds.clear();
	pos = t;
}

Vertex::Vertex(double x,double y,double z)
{
	id = -97;
	neighborVertexIds.clear();
	pos = Vec3d(x,y,z);
}


Vertex::~Vertex(void)
{
}

void Vertex::addNeighborVertex(int num){
	neighborVertexIds.insert(num);
}

void Vertex::delNeighborVertex(int num){
	neighborVertexIds.erase(num);
}

bool Vertex::hasNeighborVertex(int num){
	return (neighborVertexIds.count(num) >  0);
}


