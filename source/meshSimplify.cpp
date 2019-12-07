#include "meshSimplify.h"
MeshSimplify::MeshSimplify(void)
{
	eHeap = new EdgeHeap();
	vGroup = new VertexGroup();
	cntFace = 0;
}


MeshSimplify::~MeshSimplify(void)
{
}

void MeshSimplify::start(){
	for(int i = 0;i < cntDelFace;i += 2){ // start edge delete
		Edge e = eHeap->getMinCost();
		Vertex* v1 = &(vGroup->group[e.v1]);
		Vertex* v2 = &(vGroup->group[e.v2]);
		Vertex v0 = e.v;
		int v0_id = vGroup->addVertex(v0);
		Vertex* pV0 = &(vGroup->group[v0_id]); // set location after edge collapse

		set<int> neighborV; // pV0's neighbors 
		neighborV.clear();
		eHeap->delEdge(e); 

		for(set<int>::iterator it = v1->neighborVertexIds.begin(); it != v1->neighborVertexIds.end();it++){
			if((*it)!=v2->id){
				eHeap->delEdge( Edge((*it),v1->id));
				vGroup->group[(*it)].delNeighborVertex(v1->id);
				neighborV.insert((*it));
			}

		}

		for(set<int>::iterator it = v2->neighborVertexIds.begin(); it != v2->neighborVertexIds.end();it++){
			if((*it)!=v1->id){
				eHeap->delEdge( Edge((*it),v2->id));
				vGroup->group[(*it)].delNeighborVertex(v2->id);
				neighborV.insert((*it));
			}

		}

		for (set<int>::iterator it = neighborV.begin();it != neighborV.end(); it++) {
			vGroup->group[(*it)].addNeighborVertex(v0_id);
			vGroup->group[v0_id].addNeighborVertex(*it);
		}
		vGroup->isDeleted[v1->id] = true;
		vGroup->isDeleted[v2->id] = true;

		for (set<int>::iterator it = neighborV.begin(); it != neighborV.end(); it++) {
			Edge e((*it),v0_id);
			calVAndCostV(e);
			eHeap->addEdge(e);
		}



	}

}

void MeshSimplify::setRatio(double _ratio){
	ratio = _ratio;
}

void MeshSimplify::input(){
	int cntv=0,cntf=0;
	char s[256];
	int tmp = 0;
	while (scanf("%s", s) == 1){
		tmp++;
		switch (s[0]) {
			case '#': fgets(s, sizeof s, stdin); break;
			case 'v': {
				cntv++;
				double x, y, z;
				scanf("%lf %lf %lf", &x, &y, &z);
				vGroup -> addVertex(Vertex(x,y,z));
				break;
			}
			case 'f': {
				cntf++;
				cntFace++;
				int a, b, c;
				scanf("%d%d%d", &a, &b, &c);
				vGroup->group[a].addNeighborVertex(b);
				vGroup->group[a].addNeighborVertex(c);
				vGroup->group[b].addNeighborVertex(a);
				vGroup->group[b].addNeighborVertex(c);
				vGroup->group[c].addNeighborVertex(a);
				vGroup->group[c].addNeighborVertex(b);
				break;
			}
			default: break;
		}
	}

	for(int i = 1;i <= vGroup->cntVertex;i++){
		for(set<int>::iterator it = vGroup->group[i].neighborVertexIds.begin();
			it != vGroup->group[i].neighborVertexIds.end();it++){
			if(i < (*it))
				break;
			
			Edge e((*it),i);
			calVAndCostV(e);
			eHeap->addEdge(e);
		}
	}
	cntDelFace = (int)((1-ratio) * cntFace);
}


void MeshSimplify::output(){
	int cnt = 0;
	int cntv=0,cntf=0;
	for(int i = 1;i <= vGroup->cntVertex;i++){
		if(vGroup->isDeleted[i])
			continue;
		Vertex* v = &vGroup->group[i];
		cnt++;
		v->id = cnt;
		printf("v %lf %lf %lf\n",v->pos.x,v->pos.y,v->pos.z);
	}	
	for(int i = 1;i <= vGroup->cntVertex;i++){
		if(vGroup->isDeleted[i])
			continue;
		Vertex* v = &(vGroup->group[i]);
		for(set<int>::iterator it1 = v->neighborVertexIds.begin();it1 != v->neighborVertexIds.end();it1++){
			if(i >= (*it1))
				continue;
			for(set<int>::iterator it2 = v->neighborVertexIds.begin();it2 != v->neighborVertexIds.end();it2++){
				if((*it1) < (*it2) && (vGroup->group[(*it1)].hasNeighborVertex(*it2))){
					printf("f %d %d %d\n",v->id,vGroup->group[(*it1)].id,vGroup->group[(*it2)].id);
					cntf++;
				}	
			}
		}

	}
}

void MeshSimplify::calVAndCostV(Edge& e){
	Matrix mat = calVertexCost(e.v1) + calVertexCost(e.v2);
	e.v = calVertexPos(e,mat);
	Vector4 X(e.v.x, e.v.y, e.v.z, 1.0);
	if (vGroup->getCommonVertexNum(e.v1, e.v2) != 2) {
		e.costV = 0;
		return;
	}
	double pri = 0;
	for (int i = 0; i < 4; i++) {
		double p = 0;
		for (int j = 0; j < 4; j++)
			p += X.v[j] * mat.mat[i][j];
		pri += p * X.v[i];
	}
	e.costV = pri;
}

Matrix MeshSimplify::calVertexCost(int _id){
	Matrix ans;
	Vertex* p = &(vGroup->group[_id]);
	for(set<int>::iterator it1 = p->neighborVertexIds.begin();it1 != p->neighborVertexIds.end();it1++){
		for(set<int>::iterator it2 = p->neighborVertexIds.begin();it2 != p->neighborVertexIds.end();it2++){
			if((*it1) < (*it2) && (vGroup->group[(*it1)].hasNeighborVertex(*it2))){
				Vertex* v1 = &(vGroup->group[(*it1)]);
				Vertex* v2 = &(vGroup->group[(*it2)]);
				Vec3d n = ( (v1->pos) - (p->pos) ).getCross( (v2->pos) - (p->pos)).getUnitVectorOfThis();
				Vector4 tmp(n.x, n.y, n.z, -(p->pos.getDot(n)));
				for(int i = 0;i < 4;i++){
					for(int j = 0;j < 4;j++){
						ans.mat[i][j] += tmp.v[i] * tmp.v[j];
					}

				}

			}

		}

	}
	return ans;

}

Vec3d MeshSimplify::calVertexPos(Edge& e,Matrix m){
	Vec3d mid = (vGroup->group[e.v1].pos + vGroup->group[e.v2].pos) / 2;
	m.mat[3][0] = 0;
	m.mat[3][1] = 0;
	m.mat[3][2] = 0;
	m.mat[3][3] = 1;
	
	Vector4 Y(0,0,0,1);
	Solve* solve = new Solve(m,Y);
	Vector4 ans = solve->getAns();
	if(ans.v[3] > Config::EPS)
		return Vec3d(ans.v[0],ans.v[1],ans.v[2]);
	else
		return mid;
}






