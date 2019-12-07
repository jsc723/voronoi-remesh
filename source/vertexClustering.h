#ifndef VERTEXCLUSTER_H
#define VERTEXCLUSTER_H

#include "edgeHeap.h"
#include "vertex.h"
#include "vertexGroup.h"
#include "matrix.h"
#include "vector4.h"
#include "vector3.h"
#include "solve.h"
#include "config.h"
#include <cstdlib>
#include <cstdio>
#include <string>
#include <queue>
#include <tuple>
#include <Eigen/Dense>
using Eigen::Matrix2d;
using Eigen::Matrix3d;
using Eigen::Matrix4d;
using Eigen::Vector2d;
using Eigen::Vector3d;
using Eigen::Vector4d;
using std::vector;
using std::tuple;
using std::queue;
using std::set;


class VCCluster;
class VCFace;
class VCEdge
{
public:
	VCFace* f1, * f2;
	VCEdge() : f1(NULL), f2(NULL) {}
	VCEdge(VCFace* f1, VCFace* f2) : f1(f1), f2(f2) {}
	bool isBoundary();
	bool operator==(const VCEdge &e) const {
		return (f1 == e.f1 && f2 == e.f2) || (f2 == e.f1 && f1 == e.f2);
	}
};

class VCFace
{
public:
	int v1, v2, v3;
	Vector3d normal;
	Matrix3d K; //metric tensor
	Matrix4d E; //quadric matrix
	VCFace(int v1, int v2, int v3, VertexGroup *vg, VCCluster* c);
	void addNeighbor(VCFace* f) {
		neighbor.insert(f);
	}

	set<VCFace*> neighbor;
	double area;
	Vector3d center;
	VCCluster* cluster;
};



class VCCluster
{
private:
	void bfs(VCFace *init, set< VCFace* >* component, set<VCFace*>* unvisited = NULL);
protected:
	virtual void addItemAdditionalAct(VCFace* f) {}
	virtual void delItemAdditionalAct(VCFace* f) {}
	static const double ENG_NULL;
public:
	VCCluster(int id);
	virtual ~VCCluster(); //WARNING: delete a VCCluster will delete all of its faces
	void addItem(VCFace* f);
	void delItem(VCFace* f);
	virtual double energy();
	virtual double energyWithItem(VCFace* f);
	virtual double energyWithoutItem(VCFace* f);
	void giveItem(VCFace* f, VCCluster *c);
	virtual Vector3d center();
	size_t size() { return items.size(); }
	bool removeUnconnected(VCCluster *cNull);
	bool connected();

	int id;
	set<VCFace*> items;
	Vector3d weightedSum; //sum(p_i*y_i)
	double area;         //sum(p_i)
	bool isNull;
	Matrix4d E;
};
struct VCIQuad {
	virtual bool solveQuadMetric(const Matrix4d& E, Vector3d& result);
};
class VCClusterQuad : public VCCluster, public VCIQuad
{
public:
	VCClusterQuad(int id) : VCCluster(id) {}
	virtual double energyWithItem(VCFace* f) override;
	virtual double energyWithoutItem(VCFace* f) override;
	virtual Vector3d center() override;
	Vector3d center(const Matrix4d& E);
};
class VCClusterAniso : public VCCluster
{
protected:
	Matrix3d sumK;
	Vector3d sumKPos;
	virtual void addItemAdditionalAct(VCFace* f) override;
	virtual void delItemAdditionalAct(VCFace* f) override;
public:
	VCClusterAniso(int id) : VCCluster(id)
	{
		sumK.setZero();
		sumKPos.setZero();
	}
	virtual double energy() override;
	virtual double energyWithItem(VCFace* f) override;
	virtual double energyWithoutItem(VCFace* f) override;
	virtual Vector3d center();
	Vector3d center(const Matrix3d& sumK, const Vector3d &sumKPos);
};

class VCClusterAnisoQuad : public VCClusterAniso, public VCIQuad
{
public:
	VCClusterAnisoQuad(int id) : VCClusterAniso(id)
	{
	}
	virtual double energyWithItem(VCFace* f) override;
	virtual double energyWithoutItem(VCFace* f) override;
	virtual Vector3d center() override;
	Vector3d center(const Matrix4d& E, const Matrix3d& sumK, const Vector3d& sumKPos);
};

struct VCPair {
	int a, b;
	VCPair(int a, int b) :
		a(a), b(b) {}
};

struct VCEdgeItemCmp {
	bool operator()(const VCPair& lhs, const VCPair& rhs) const {
		if (lhs.a < rhs.a)
			return true;
		if (lhs.a > rhs.a)
			return false;
		return lhs.b < rhs.b;
	}
	bool operator()(const VCEdge& lhs, const VCEdge& rhs) const {
		if (lhs.f1 < rhs.f1)
			return true;
		if (lhs.f1 > rhs.f1)
			return false;
		return lhs.f2 < rhs.f2;
	}
};


class VertexClustering
{
public:
	struct Options {
		size_t numCluster;
		bool adaptive;
		bool anisotropic;
		bool quadMetric;
		bool validation;
		int ringLevel;
		int seed;
	};
	VertexClustering(Options &opts);
	~VertexClustering(void);
	void input(const string& filePath);
	void start();
	void build();
	void output(const string& filePath);
	void stat();
private:
	VertexGroup* vGroup;
	VertexGroup* vGroupNew;
	vector<VCCluster *> clusters;
	map<VCPair, VCEdge, VCEdgeItemCmp> edgeMap;
	queue<VCEdge> boundary;
	size_t cntFace;
	double angle_min;
	double angle_less_30;
	Options opts;
	
	void initClusters(vector<VCEdge>& edges);
	void pushNeighbors(VCFace* f, VCEdge& e, queue<VCEdge> *q);
	double totalEnergy();

	void calcCurvature(vector<VCFace *> &faces);
	set<int> collectRings(VCFace* f);
	vector<double> leastSquareFitting(const vector<Vector3d> &data);
	double localCurvature(const vector<double>& poly, Vector2d& egs, Matrix2d& evs);

	double quadricCost(const Matrix4d& E, const Vector3d& v) {
		Vector4d v4;
		v4 << v(0), v(1), v(2), 1;
		return v4.transpose() * E * v4;
	}
	void localRefine(const set<int> &cset);
	deque<int> sortNeighbor(const set<int>& neighbor);
	void computeClusters(bool constrain);
	bool removeUnconnected();
	bool transferLeadsToDisconnect(VCCluster* c1, VCFace* f);
	bool localEnergyRelease(VCCluster* c1, VCCluster* c2, VCEdge& e, queue<VCEdge>* q2, bool constrain);
};


#endif