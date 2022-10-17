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
	VCFrontEnd::exec(argc, argv);
	return 0;
}
