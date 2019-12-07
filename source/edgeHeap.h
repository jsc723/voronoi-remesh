#ifndef EDGEHEAP_H
#define EDGEHEAP_H
#include "config.h"
#include "edge.h"
#include <queue>
#include <vector>
#include <map>
#include <iostream>
using namespace std;
class EdgeHeap
{
public:
	EdgeHeap(void);
	~EdgeHeap(void);
	struct cmp{
		bool operator() (Edge X, Edge Y){
			return X.costV > Y.costV;
		}
	};
	std::priority_queue<Edge,std::vector<Edge>,cmp> pq; // used to store the edges, ordered by delta from small to large
	map<pair<int, int>, int> mapEdgeToID; // create the mapping from edge to vertex
	bool isDeleted[Config::MAX_NUM_EDGES+1]; // keep track of which edges are deleted
	int cntEdge; // number of edges
	void addEdge(Edge&); 
	void delEdge(Edge);
	Edge getMinCost(); // delete the edge with smallest delta
};
#endif
