#include "header.h"
#include "meshSimplify.h"
#include "config.h"
#include "matrix.h"
#include "vector4.h"
#include "solve.h"
#include <iostream>
#include <time.h>
#include "vertexClustering.h"
#include "VCFrontEnd.h"

using namespace std;

int main(int argc,char* argv[]){	
	if(argc < 2){
		cout << "error input " << endl;
		return -1;
	}else if(string(argv[1]) == "voronoi"){
		VCFrontEnd::exec(argc - 1, &argv[1]);
	}else if(string(argv[1]) == "meshInterpolate"){
		cout << "not supported in this open source version" << endl;
	}



	return 0;
}
