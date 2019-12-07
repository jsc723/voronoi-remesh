#ifndef MESHSIMPLIFY_H
#define MESHSIMPLIFY_H
#include "edgeHeap.h"
#include "vertexGroup.h"
#include "matrix.h"
#include "vector4.h"
#include "solve.h"
#include "config.h"
#include <cstdlib>
#include <cstdio>
#include <string>
using namespace std;
class MeshSimplify
{
	double ratio;
	int cntFace,cntDelFace;
	EdgeHeap* eHeap;
	VertexGroup* vGroup;
public:
	MeshSimplify(void);
	~MeshSimplify(void);

	void setRatio(double);
	void setLeftFaceNum(int);

	void input();
	void start();
	void output();

	Matrix calVertexCost(int);
	Vec3d calVertexPos(Edge&,Matrix);
	
	void calVAndCostV(Edge&);
};

#endif