#include "vertexClustering.h"
#include "vector3.h"
#include <set>
#include <ctime>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <limits>

using namespace Eigen;
using namespace std;

/* some helper functions */

static Vec3d VCEigenToMyvec(const Vector3d& v) {
	return Vec3d(v.x(), v.y(), v.z());
}
static Vector3d VCMyvecToEigen(const Vec3d& v) {
	return Vector3d(v.x, v.y, v.z);
}

static Vector3d operator- (const Vec3d& v1, const Vector3d& v2) {
	return Vector3d(v1.x - v2.x(), v1.y - v2.y(), v1.z - v2.z());
}
static Vector3d operator- (const Vector3d & v1, const Vec3d & v2) {
	return Vector3d(v1.x() - v2.x, v1.y() - v2.y, v1.z() - v2.z);
}
static double angleBetween(Vec3d v1, Vec3d v2) {
	const double PI = atan(1) * 4;
	Vec3d u1(v1.getUnitVectorOfThis()), u2(v2.getUnitVectorOfThis());
	return acos(u1.getDot(u2))* 180.0 / PI; //in degrees
}


/* ----------------VCEdge ----------------*/

bool VCEdge::isBoundary() {
	if (!f2 || !f1) return false;
	return f1->cluster != f2->cluster;
}

/*----------------- VCFace ------------------*/

VCFace::VCFace(int v1, int v2, int v3, VertexGroup* vg, VCCluster* c)
	:v1(v1), v2(v2), v3(v3), cluster(c)
{
	Vertex* p1 = &vg->group[v1];
	Vertex* p2 = &vg->group[v2];
	Vertex* p3 = &vg->group[v3];
	Vector3d u = VCMyvecToEigen(p2->pos - p1->pos);
	Vector3d v = VCMyvecToEigen(p3->pos - p1->pos);
	Vector3d n(u.cross(v)), a(VCMyvecToEigen(p1->pos));
	normal = n.normalized();
	area = 0.5 * sqrt(n.dot(n));
	center = VCMyvecToEigen((p1->pos + p2->pos + p3->pos) / 3);
	n.normalize();
	Vector4d plane;
	plane << n.x() , n.y() , n.z() ,-n.dot(a); 
	E = plane * plane.transpose();
	c->addItem(this);
}


/*---------------- VCCluster ---------------*/

VCCluster::VCCluster(int id):id(id), area(0), isNull(false)
{
	E.setZero();
}
VCCluster::~VCCluster() {
	for (VCFace *f : items)
		delete f;
}
/* bfs helper function for removeUnconnected() and connected() */
void VCCluster::bfs(VCFace* init, set<VCFace *> *component, set<VCFace*> *unvisited)
{
	queue<VCFace*> q;
	q.push(init);
	while (!q.empty()) {
		VCFace *f = q.front();
		q.pop();
		if (component->find(f) == component->end()) {
			if(unvisited)
				unvisited->erase(f);
			component->insert(f);
			for (VCFace* n : f->neighbor) {
				if (n->cluster == this)
					q.push(n);
			}
		}
	}
}
/* 
	Find all connected components in the cluster. If there are more than one, 
	then transfer all components except the largest one to cNull.
*/
bool VCCluster::removeUnconnected(VCCluster* cNull) {
	vector< set<VCFace*> * > components;
	set<VCFace*> unvisited(items);
	while (!unvisited.empty()) {
		set<VCFace*>* comp = new set<VCFace*>();
		bfs(*unvisited.begin(), comp, &unvisited);
		components.push_back(comp);
	}
	if (components.size() <= 1) {
		if(components.size() == 1)
			delete components[0];
		return false;
	}
	set<VCFace*>* largest = *max_element(components.begin(), components.end(),
		[](set<VCFace*> * s1, set<VCFace*> * s2) -> bool { return s1->size() < s2->size(); });
	for (set<VCFace*>* s : components) {
		if (s != largest) {
			for (VCFace* item : *s) {
				giveItem(item, cNull);
			}
			delete s;
		}
	}
	delete largest;
	return true;
}

const double VCCluster::ENG_NULL = 0;

void VCCluster::addItem(VCFace* f) {
	if (items.find(f) == items.end()) {
		items.insert(f);
		area += f->area;
		weightedSum = weightedSum + f->center * f->area;
		f->cluster = this;
		E += f->E;
		addItemAdditionalAct(f);
	}
}
void VCCluster::delItem(VCFace* f) {
	auto it = items.find(f);
	if (it != items.end()) {
		items.erase(it);
		area -= f->area;
		weightedSum = weightedSum - f->center * f->area;
		f->cluster = NULL;
		E -= f->E;
		delItemAdditionalAct(f);
	}
}

void VCCluster::giveItem(VCFace* f, VCCluster *c) {
	delItem(f);
	c->addItem(f);
}
Vector3d VCCluster::center() {
	if (area == 0)
		return Vector3d(0,0,0);
	return weightedSum / area;
}
double VCCluster::energy() {
	Vector3d zi = center();
	return area == 0 ? ENG_NULL : area * zi.dot(zi) - 2 * zi.dot(weightedSum);
}
double VCCluster::energyWithItem(VCFace* f) {
	Vector3d ws = weightedSum + f->center * f->area;
	double a = area + f->area;
	Vector3d zi = center();
	return a * zi.dot(zi) - 2 * zi.dot(ws);
}
double VCCluster::energyWithoutItem(VCFace* f) {
	Vector3d ws = weightedSum - f->center * f->area;
	double a = area - f->area;
	Vector3d zi = center();
	return a == 0 ? ENG_NULL : a * zi.dot(zi) - 2 * zi.dot(ws);
}

bool VCCluster::connected() {
	if (items.size() < 2) return true;
	set<VCFace*> c;
	bfs(*items.begin(), &c);
	return c.size() == items.size();
}

/* VCIQuad */

bool VCIQuad::solveQuadMetric(const Matrix4d& E, Vector3d& result) {
	Matrix4d E1;
	E1.row(0) = E.row(0);
	E1.row(1) << E(0, 1), E(1, 1), E(1, 2), E(1, 3);
	E1.row(2) << E(0, 2), E(1, 2), E(2, 2), E(2, 3);
	E1.row(3) << 0, 0, 0, 1;
	Vector4d b, x;
	b << 0, 0, 0, 1;
	x = E1.colPivHouseholderQr().solve(b);
	double relative_error = (E1 * x - b).norm() / b.norm();
	result(0) = x(0);
	result(1) = x(1);
	result(2) = x(2);
	return relative_error < 1e-6;
}

/* VCClusterQuad */

double VCClusterQuad::energyWithItem(VCFace* f) {
	Matrix4d Eo = E + f->E;
	Vector3d ws = weightedSum + f->center * f->area;
	double a = area + f->area;
	Vector3d zi = center(Eo);
	return a * zi.dot(zi) - 2 * zi.dot(ws);
}
double VCClusterQuad::energyWithoutItem(VCFace* f) {
	if (items.size() <= 1)
		return ENG_NULL;
	Matrix4d Eo = E - f->E;
	Vector3d ws = weightedSum - f->center * f->area;
	double a = area - f->area;
	Vector3d zi = center(Eo);
	return a == 0 ? ENG_NULL : a * zi.dot(zi) - 2 * zi.dot(ws);
}

Vector3d VCClusterQuad::center() {
	return center(E);
}
Vector3d VCClusterQuad::center(const Matrix4d& E) {
	Vector3d x;
	if (solveQuadMetric(E, x))
		return x;
	return VCCluster::center();
}

/* VCClusterAniso */

Vector3d VCClusterAniso::center() {
	return center(sumK, sumKPos);
}
Vector3d VCClusterAniso::center(const Matrix3d& sumK, const Vector3d &sumKPos)
{
	bool valid;
	Matrix3d sumKInv;
	double det;
	sumK.computeInverseAndDetWithCheck(sumKInv, det, valid);
	if (!valid) {
		return Vector3d(0, 0, 0);
	}
	return sumKInv * sumKPos;
}
void VCClusterAniso::addItemAdditionalAct(VCFace* f) {
	sumK = sumK + f->K;
	sumKPos = sumKPos + f->K * f->center;
}
void VCClusterAniso::delItemAdditionalAct(VCFace* f) {
	sumK = sumK - f->K;
	sumKPos = sumKPos - f->K * f->center;
}

double VCClusterAniso::energy() {
	Vector3d zi = center();
	return zi.dot(sumK * zi) - 2 * zi.dot(sumKPos);
}
double VCClusterAniso::energyWithItem(VCFace* f) {
	Matrix3d sK = sumK + f->K;
	Vector3d sKP = sumKPos + f->K * f->center;
	Vector3d zi = center(sK, sKP);
	return zi.dot(sK * zi) - 2 * zi.dot(sKP);
}
double VCClusterAniso::energyWithoutItem(VCFace* f) {
	if (items.size() <= 1)
		return ENG_NULL;
	Matrix3d sK = sumK - f->K;
	Vector3d sKP = sumKPos - f->K * f->center;
	Vector3d zi = center(sK, sKP);
	return zi.dot(sK * zi) - 2 * zi.dot(sKP);
}

/* VCClusterAnisoQuad */

double VCClusterAnisoQuad::energyWithItem(VCFace* f) {
	Matrix3d sK = sumK + f->K;
	Vector3d sKP = sumKPos + f->K * f->center;
	Matrix4d Eo = E + f->E;
	Vector3d zi = center(Eo, sK, sKP);
	return zi.dot(sK * zi) - 2 * zi.dot(sKP);
}
double VCClusterAnisoQuad::energyWithoutItem(VCFace* f) {
	if (items.size() <= 1)
		return ENG_NULL;
	Matrix3d sK = sumK - f->K;
	Vector3d sKP = sumKPos - f->K * f->center;
	Matrix4d Eo = E - f->E;
	Vector3d zi = center(Eo, sK, sKP);
	return zi.dot(sK * zi) - 2 * zi.dot(sKP);
}

Vector3d VCClusterAnisoQuad::center() {
	return center(E, sumK, sumKPos);
}
Vector3d VCClusterAnisoQuad::center(const Matrix4d& E, const Matrix3d& sumK, const Vector3d& sumKPos) {
	Vector3d x;
	if (solveQuadMetric(E, x))
		return x;
	return VCClusterAniso::center(sumK, sumKPos);
}

/* ----------------- VertexClustering ------------------*/

VertexClustering::VertexClustering(VertexClustering::Options &opts): cntFace(0)
{
	srand(opts.seed);
	if (opts.anisotropic) {
		opts.adaptive = true;
	}
	this->opts = opts;
	vGroup = new VertexGroup();
	vGroupNew = NULL;
	angle_min = 180;
	angle_less_30 = 100;
	clusters.resize(opts.numCluster + 1); 
	if (opts.anisotropic) {
		if (opts.quadMetric) {
			for (int i = 0; i < clusters.size(); i++)
				clusters[i] = new VCClusterAnisoQuad(i);
		}
		else {
			for (int i = 0; i < clusters.size(); i++)
				clusters[i] = new VCClusterAniso(i);
		}
	}
	else
		if (opts.quadMetric) {
			for (int i = 0; i < clusters.size(); i++)
				clusters[i] = new VCClusterQuad(i);
		}
		else {
			for (int i = 0; i < clusters.size(); i++)
				clusters[i] = new VCCluster(i);
		}
	clusters[0]->isNull = true;
}

VertexClustering::~VertexClustering(void)
{
	for (int i = 0; i < clusters.size(); i++)
		delete clusters[i];
	delete vGroup;
	if (vGroupNew)
		delete vGroupNew;
}

/* sort three number in ascending order */
template< typename T >
inline void VCSort3(T& a, T& b, T& c) {
	if (b > c) std::swap(b, c);
	if (a > b) std::swap(a, b);
	if (b > c) std::swap(b, c);
}

/* Read the mesh from the file */
void VertexClustering::input(const string& filePath)
{
	FILE* f = fopen(filePath.c_str(), "r");
	int cntv = 0, cntf = 0;
	char s[256], _;
	vector<VCEdge> edges;
	vector<VCFace *> faces;

	while (fgets(s, sizeof s, f) != NULL) {
		switch (s[0]) {
		case '#': break; 
		case 'v': {
			cntv++;
			double x, y, z;
			sscanf(s, "%c%lf%lf%lf", &_, &x, &y, &z);
			vGroup->addVertex(Vertex(x, y, z));
			break;
		}
		case 'f': {
			cntf++;
			cntFace++;
			int a, b, c;
			sscanf(s, "%c%d%d%d",&_, &a, &b, &c);
			VCSort3(a, b, c);
			vGroup->group[a].addNeighborVertex(b);
			vGroup->group[a].addNeighborVertex(c);
			vGroup->group[b].addNeighborVertex(a);
			vGroup->group[b].addNeighborVertex(c);
			vGroup->group[c].addNeighborVertex(a);
			vGroup->group[c].addNeighborVertex(b);
			VCFace* f = new VCFace(a, b, c, vGroup, clusters[0]); //face pointers tracked by clusters
			faces.push_back(f);
			VCPair pairs[3] = { VCPair(a,b), VCPair(a,c), VCPair(b,c) };
			for (int i = 0; i < 3; ++i) {
				auto iter = edgeMap.find(pairs[i]);
				if (iter != edgeMap.end()) {
					iter->second.f2 = f;
					iter->second.f1->addNeighbor(f);
					f->addNeighbor(iter->second.f1);
					edges.push_back(iter->second);
				}
				else {
					edgeMap[pairs[i]] = VCEdge(f, NULL);
				}
			}
			break;
		}
		default: break;
		}
	}
	if(opts.adaptive)
		calcCurvature(faces); //only uniform remeshing can skip this step
	initClusters(edges);
	for (auto e = edges.begin(); e != edges.end(); ++e) {
		if (e->isBoundary())
			boundary.push(*e);
	}
}
/* 
	Initialize clusters so that each cluster is associate to at least one triangle.
*/
void VertexClustering::initClusters(vector<VCEdge>& edges) {
	int n = edges.size();
	set<VCFace*> faces;
	if (opts.numCluster > n)
		opts.numCluster = n;
	while (faces.size() < opts.numCluster) {
		while (true) {
			int r = rand() % n;
			VCFace* f = edges[r].f1;
			if (faces.find(f) == faces.end()) {
				faces.insert(f);
				break;
			}
		}
	}
	vector<VCFace*> faces_vec(faces.begin(), faces.end());
	for (int i = 0; i < opts.numCluster; ++i) {
		clusters[0]->giveItem(faces_vec[i], clusters[i + 1]);
	}
	
}

/* push all of the edges of f except e in to q */
void VertexClustering::pushNeighbors(VCFace *f, VCEdge& e, queue<VCEdge>* q)
{
	VCEdge e1(edgeMap[VCPair(f->v1, f->v2)]);
	VCEdge e2(edgeMap[VCPair(f->v1, f->v3)]);
	VCEdge e3(edgeMap[VCPair(f->v2, f->v3)]);
	if (e1 == e) {
		q->push(e2);
		q->push(e3);
	}
	else if (e2 == e) {
		q->push(e1);
		q->push(e3);
	}
	else {
		q->push(e1);
		q->push(e2);
	}
}
/* Return total energy of all clusters */
double VertexClustering::totalEnergy() {
	double energy = 0;
	for (int i = 0; i < clusters.size(); ++i) {
		energy += clusters[i]->energy();
	}
	return energy;
}

/* Determine whether delete f from c1 will make c1 to be disconnected */
bool VertexClustering::transferLeadsToDisconnect(VCCluster* c1, VCFace* f) {
	c1->delItem(f);
	bool ret = c1->connected();
	c1->addItem(f);
	return !ret;
}
/* 
	redistribute the 2 triangles which shares the edge e 
	according to the minimum energy.
*/
bool VertexClustering::localEnergyRelease(VCCluster* c1, VCCluster* c2, VCEdge& e, queue<VCEdge> *q2, bool constrain) {
	double L0, L1, L2;
	if (c1->isNull) goto localEnergyRelease_c1_c2;
	if (c2->isNull) goto localEnergyRelease_c2_c1;
	L0 = c1->energy() + c2->energy();
	L1 = c1->energyWithItem(e.f2) + c2->energyWithoutItem(e.f2);
	L2 = c1->energyWithoutItem(e.f1) + c2->energyWithItem(e.f1);
	double minEnergy;
	if (constrain) {
		bool c1_ok = transferLeadsToDisconnect(c1, e.f1);
		bool c2_ok = transferLeadsToDisconnect(c2, e.f2);
		if (!c2_ok) 
			L1 = numeric_limits<double>::max();
		if (!c1_ok) 
			L2 = numeric_limits<double>::max();
	}
	minEnergy = min(min(L0, L1), L2);
	if (minEnergy == L1) 
		goto localEnergyRelease_c2_c1;
	else if (minEnergy == L2) 
		goto localEnergyRelease_c1_c2;
	q2->push(e);
	return false;

localEnergyRelease_c2_c1:
	c2->giveItem(e.f2, c1);
	pushNeighbors(e.f2, e, q2);
	return true;
localEnergyRelease_c1_c2:
	c1->giveItem(e.f1, c2);
	pushNeighbors(e.f1, e, q2);
	return true;
}

/*
	Run the vertex clustering algorithm
	Should initialize clusters and boundary edges before calling this function.
*/
void VertexClustering::computeClusters(bool constrain)
{
	queue<VCEdge> boundary_next;
	queue<VCEdge> *q1 = &boundary, *q2 = &boundary_next;
	int modification;
	int iteration = 0;
	clock_t start = clock();
	printf("---start clustering, with total energy: %.12lf---\n", totalEnergy());
	do {
		modification = 0;
		set<VCEdge, VCEdgeItemCmp> tested;
		while (!q1->empty()) {
			VCEdge e = q1->front();
			q1->pop();
			if (e.isBoundary() && tested.find(e) == tested.end()) {
				tested.insert(e);
				VCCluster *c1 = e.f1->cluster;
				VCCluster *c2 = e.f2->cluster;
				if (localEnergyRelease(c1, c2, e, q2, constrain)) {
					modification++;
				}
			}
		}
		printf("iteration: %d, modification: %d, ", ++iteration, modification);
		printf("total energy: %.12lf\n", totalEnergy());
		std::swap(q1, q2);
	} while (modification > 0);
	if (boundary.empty())
		boundary = boundary_next;

	printf("Done, uses %lf s for clustering.\n", (double)(clock() - start) / CLOCKS_PER_SEC);
	printf("---total energy: %.12lf---\n", totalEnergy());
}

bool VertexClustering::removeUnconnected() {
	bool notDone = false;
	for (int i = 1; i < clusters.size(); i++) {
		if (clusters[i]->removeUnconnected(clusters[0])) {
			notDone = true;
		}
	}
	return notDone;
}
void VertexClustering::start() {
	printf("--- first clustering iteration ---\n");
	computeClusters(false); 
	if (opts.validation) {
		if (removeUnconnected()) {
			printf("---second clustering iteration---\n");
			computeClusters(true);
		}
	}
}

/*
	Rebuild the mesh from the clusters 
*/
void VertexClustering::build()
{
	printf("start to rebuild the mesh from the clusters\n");
	clock_t start = clock();
	if (vGroupNew)
		delete vGroupNew;
	vGroupNew = new VertexGroup();
	for (int i = 1; i < clusters.size(); i++) {
		vGroupNew->addVertex(Vertex(VCEigenToMyvec(clusters[i]->center())));
		if (clusters[i]->size() == 0)
			vGroupNew->delVertex(i);
	}

	while(boundary.size() > 0) {
		VCEdge e = boundary.front();
		boundary.pop();
		if (e.isBoundary()) {
			int v1 = e.f1->cluster->id;
			int v2 = e.f2->cluster->id;
			vGroupNew->group[v1].addNeighborVertex(v2);
			vGroupNew->group[v2].addNeighborVertex(v1);
		}
	}
	//refine
	for (int i = 1; i <= vGroup->cntVertex; i++) {
		if (vGroup->isDeleted[i])
			continue;
		set<int> cset;
		Vertex* v = &vGroup->group[i];
		for (int u : v->neighborVertexIds) {
			VCPair pair(min(i, u), max(i, u));
			auto edgeIt = edgeMap.find(pair);
			if (edgeIt != edgeMap.end()) {
				if(edgeIt->second.f1)
					cset.insert(edgeIt->second.f1->cluster->id);
				if(edgeIt->second.f2)
					cset.insert(edgeIt->second.f2->cluster->id);
			}
		}
		if (cset.size() >= 4) {
			localRefine(cset);
		}
	}
	printf("Done, uses %lf s for building.\n", (double)(clock() - start) / CLOCKS_PER_SEC);
}
/* 
	given the vertices of a polygon, add edges between them so that the polygon is divided into triangles 
	*/
void VertexClustering::localRefine(const set<int> &cset)
{
	//build a min heap
	deque<int> ring = sortNeighbor(cset);
	if (ring.size() != cset.size())
		return; //not a maniford, cannot refine
	typedef tuple<int, int, double> CostPair;
	auto cmp = [](const CostPair e1, const CostPair & e2) -> bool {return get<2>(e1) < get<2>(e2); };
	set<CostPair, decltype(cmp)> candidates(cmp);
	Matrix4d E;
	E.setZero();
	for (int c : cset)
		E += clusters[c]->E;
	for (int a : cset) {
		for (int b : cset) {
			if (a < b && !vGroupNew->group[a].hasNeighborVertex(b)) {
				Vec3d v1 = vGroupNew->group[a].pos;
				Vec3d v2 = vGroupNew->group[b].pos;
				Vector3d m = VCMyvecToEigen(v1 + v2) / 2;
				double cost = quadricCost(E, m);
				candidates.insert(make_tuple(a, b, cost));
			}
		}
	}
	int r = ring.size();
	if (candidates.size() <= r - (int)3)
		return;
	auto find_order = [&ring](int i) -> int
		{ return find(ring.begin(), ring.end(), i) - ring.begin(); };
	while (!candidates.empty()) {
		auto first = candidates.begin();
		CostPair pair = *first;
		candidates.erase(first);
		int a = get<0>(pair), b = get<1>(pair);
		vGroupNew->group[a].addNeighborVertex(b);
		vGroupNew->group[b].addNeighborVertex(a);
		int ai = find_order(a), bi = find_order(b);
		if (ai > bi)
			swap(ai, bi);
		auto iter = candidates.begin();
		while (iter != candidates.end()) {
			CostPair e = *iter;
			int ci = find_order(get<0>(e)), di = find_order(get<1>(e));
			if (ci > di) {
				swap(ci, di);
			}
			if ((ci < ai && di < bi) || (ci > ai && di > bi)) {
				iter = candidates.erase(iter);
			}
			else {
				++iter;
			}
		}
	}
}

/*
	given the vertices of a polygon, sort them and return a list 
	such that each pair of list[i] and list[i+1] share an edge
*/
deque<int> VertexClustering::sortNeighbor(const set<int>& ring) {
	set<int> unsorted(ring);
	int first = *(unsorted.begin());
	unsorted.erase(first);
	deque<int> q;
	q.push_back(first);
	while (!unsorted.empty()) {
		bool modified = false;
		for (int n : unsorted) {
			if (vGroupNew->group[q.front()].hasNeighborVertex(n)) {
				q.push_front(n);
				unsorted.erase(n);
				modified = true;
				break;
			}
			if (vGroupNew->group[q.back()].hasNeighborVertex(n)) {
				q.push_back(n);
				unsorted.erase(n);
				modified = true;
				break;
			}
		}
		if (!modified)
			break;
	}
	return q;
}
/*
	Write the new mesh into the file
*/
void VertexClustering::output(const string &filePath)
{
	FILE *f = fopen(filePath.c_str(), "w");
	if (!vGroupNew) {
		printf("please run build() before output()\n");
		return;
	}
	int cnt = 0;
	int cntf = 0;
	int cnt_less_30 = 0;
	for (int i = 1; i <= vGroupNew->cntVertex; i++) {
		if (vGroupNew->isDeleted[i])
			continue;
		Vertex* v = &vGroupNew->group[i];
		cnt++;
		v->id = cnt;
		fprintf(f, "v %lf %lf %lf\n", v->pos.x, v->pos.y, v->pos.z);
	}
	for (int i = 1; i <= vGroupNew->cntVertex; i++) {
		if (vGroupNew->isDeleted[i])
			continue;
		Vertex* v = &(vGroupNew->group[i]);
		for (set<int>::iterator it1 = v->neighborVertexIds.begin(); it1 != v->neighborVertexIds.end(); it1++) {
			if (i >= (*it1))
				continue;
			for (set<int>::iterator it2 = v->neighborVertexIds.begin(); it2 != v->neighborVertexIds.end(); it2++) {
				if ((*it1) < (*it2) && (vGroupNew->group[*it1].hasNeighborVertex(*it2))) {
					Vertex* v1 = &vGroupNew->group[*it1];
					Vertex* v2 = &vGroupNew->group[*it2];
					double a1 = angleBetween(v1->pos - v->pos, v2->pos - v->pos);
					double a2 = angleBetween(v->pos - v1->pos, v2->pos - v1->pos);
					double a3 = angleBetween(v->pos - v2->pos, v1->pos - v2->pos);
					VCSort3(a1, a2, a3);
					if (a1 < angle_min) angle_min = a1;
					if (a1 < 30) cnt_less_30++;
					if (a2 < 30) cnt_less_30++;
					fprintf(f, "f %d %d %d\n", v->id, vGroupNew->group[(*it1)].id, vGroupNew->group[(*it2)].id);
					cntf++;
				}
			}
		}
	}
	angle_less_30 = 100.0 * cnt_less_30 / (3.0 * cntf);
}

void VertexClustering::stat() {
	printf("angle min: %lf\n", angle_min);
	printf("angle less than 30: %lf%%\n", angle_less_30);
}

/* ------------- Adaptive Remeshing -------------- */
/*
	Given a face (triangle) f, find all vertices that is reachable 
	in opts.ringLevel (default 2) steps from the triangle,
	including the vertices of itself.
*/
set<int> VertexClustering::collectRings(VCFace* f)
{
	set<int> rings;
	queue< tuple<int, int> > q;
	q.push(make_tuple(f->v1, 0));
	q.push(make_tuple(f->v2, 0));
	q.push(make_tuple(f->v3, 0));
	while (!q.empty()) {
		tuple<int, int> t = q.front();
		q.pop();
		int d = get<1>(t);
		if (d > opts.ringLevel)
			break;
		int v = get<0>(t);
		if (rings.find(v) == rings.end()) {
			rings.insert(v);
			for (int neighbor : vGroup->group[v].neighborVertexIds) {
				q.push(make_tuple(neighbor, d + 1));
			}
		}
	}
	return rings;
}
/* From CS184 Project 3, given a normal vector, build a orthogonal coordinate space
	 where z direction is the normal */
void make_coord_space(Matrix3d &o2w, const Vector3d & n) {
	Vector3d z = Vector3d(n);
	Vector3d h = z;
	if (fabs(h.x()) <= fabs(h.y()) && fabs(h.x()) <= fabs(h.z())) h.x() = 1.0;
	else if (fabs(h.y()) <= fabs(h.x()) && fabs(h.y()) <= fabs(h.z())) h.y() = 1.0;
	else h.z() = 1.0;
	
	z = z.normalized();
	Vector3d y = h.cross(z);
	y = y.normalized();
	Vector3d x = z.cross(y);
	x = x.normalized();
	
	o2w.col(0) = x;
	o2w.col(1) = y;
	o2w.col(2) = z;
}

/*
	Using Least Square to find a 2d polynormial fit of the given data. 
	The model is: w1 * x + w2 * y + w3 * x * x + w4 * x * y + w5 * y * y  (no constant term)
*/
vector<double> VertexClustering::leastSquareFitting(const vector<Vector3d>& data)
{
	const int SIZE = 5;
	MatrixXd A(data.size(), SIZE);
	VectorXd w(SIZE), b(data.size());
	for (int i = 0; i < data.size(); i++) {
		A.row(i) << data[i].x(),
			data[i].y(),
			data[i].x() * data[i].x(),
			data[i].x()* data[i].y(),
			data[i].y()* data[i].y();
		b(i) = data[i].z();
	}
	MatrixXd At(A.transpose());
	w = (At * A).inverse() * At * b;
	vector<double> ret(SIZE);
	for (int i = 0; i < ret.size(); i++)
		ret[i] = w(i);
	return ret;
}

/*
	Given the polynomial fitting of the surface around a face, calculate
	the curvature of the polynomial.
*/
double VertexClustering::localCurvature(const vector<double> & poly, Vector2d &egs, Matrix2d &evs)
{
	double sa1 = poly[0] * poly[0];
	double sa2 = poly[1] * poly[1];
	double b = 1 / sqrt(sa1 + sa2 + 1);
	double E = 1 + sa1;
	double F = poly[1] * poly[0];
	double G = 1 + sa2;
	double e = 2 * poly[2] * b;
	double f = poly[3] * b;
	double g = 2 * poly[4] * b;
	Matrix2d m1, m2;
	m1 << -e, -f, -f, -g;
	m2 << E, F, F, G;
	Matrix2d m = m1 * m2.inverse();
	Eigen::EigenSolver<Matrix2d> solver;
	solver.compute(m, true);
	egs = solver.eigenvalues().real();
	evs = solver.eigenvectors().real();
	return egs.dot(egs);
}

/*
	Compute curvature and Riemannian metric for each triangle
*/
void VertexClustering::calcCurvature(vector<VCFace*> & faces)
{
	printf("start to calculate curvature for %lu faces...\n", faces.size());
	clock_t start = clock();
	int cnt = 0, n10 = faces.size() / 10, percent = 1;
	for (VCFace* f : faces) {
		set<int> rings(collectRings(f));
		vector<Vector3d> pts(rings.size());
		int i = 0;
		Vector3d origin = f->center;
		for (int id : rings)
			pts[i++] = vGroup->group[id].pos - origin;
		Matrix3d o2w;
		make_coord_space(o2w, f->normal);
		Matrix3d w2o = o2w.transpose();
		for (i = 0; i < pts.size(); i++) {
			pts[i] = w2o * pts[i]; //now pts are in object coordinate
		}
		vector<double> poly(leastSquareFitting(pts));
		Matrix2d eVecs;
		Vector2d eValues;
		double weight = localCurvature(poly, eValues, eVecs);
		f->area *= weight;       //update density based on curvature
		if (opts.anisotropic) {  
			//D1, D2 are the priciple directions in xyz coordinate, their eigenvalues are principle curvatures
			Vector3d D1 = (o2w.col(0) * eVecs(0, 0) + o2w.col(1) * eVecs(1, 0)).normalized();
			Vector3d D2 = (o2w.col(0) * eVecs(0, 1) + o2w.col(1) * eVecs(1, 1)).normalized();
			Matrix3d Rotation;
			Rotation.col(0) = D1;
			Rotation.col(1) = D2;
			Rotation.col(2) = D1.cross(D2);
			Matrix3d Stretch;
			Stretch << fabs(eValues(0)), 0, 0,
				0, fabs(eValues(1)), 0,
				0, 0, max(fabs(eValues(0)), fabs(eValues(1)));
			f->K = Rotation * Stretch * Rotation.transpose(); //compute the Riemannian metric
		}
		if (++cnt > n10) {
			printf("%%%d...\n", 10 * percent);
			n10 = faces.size() / 10.0 * ++percent;
		}
	}
	printf("Done, uses %lf s for curvature.\n", (double)(clock() - start) / CLOCKS_PER_SEC);
}

