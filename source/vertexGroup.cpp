#include "vertexGroup.h"
VertexGroup::VertexGroup(void)
{
	cntVertex = 0;
	isDeleted = new bool[Config::MAX_NUM_VERTICES];
	for(int i = 0;i < Config::MAX_NUM_VERTICES;i++)
		isDeleted[i] = false;
}
VertexGroup::~VertexGroup(void)
{
}

int VertexGroup::addVertex(Vertex p){
	cntVertex++;
	p.id = cntVertex;
	group[cntVertex] = p;
	return cntVertex;
}

void VertexGroup::delVertex(int _id){
	if(_id >= Config::MAX_NUM_VERTICES){
		return;
	}	
	isDeleted[_id] = true;

	for(set<int>::iterator it = group[_id].neighborVertexIds.begin();it != group[_id].neighborVertexIds.end();it++){
		group[(*it)].delNeighborVertex(_id);
	}

}

int VertexGroup::getCommonVertexNum(int u,int v){
	int cnt = 0;
	for (set<int>::iterator it = group[u].neighborVertexIds.begin();
		it != group[u].neighborVertexIds.end();it++){
			if(group[v].hasNeighborVertex(*it)){
				cnt++;
			}

	}
	return cnt;
}
