#if defined(_WIN32)
#include <Windows.h>
#elif defined(__unix__)
#include <unistd.h>
#else
#error Platform not supported
#endif
#include <time.h>
#include <sstream>
#include <algorithm>
#include "VCFrontEnd.h"
#include "meshSimplify.h"

using namespace std;

int VCFrontEnd::exec(int argc, char** argv) {
	VCTools::InputParser input(argc, argv);
	VertexClustering::Options opts;
	opts.adaptive = false;
	opts.anisotropic = false;
	opts.numCluster = 400;
	opts.quadMetric = false;
	opts.validation = false;
	opts.ringLevel = 2;
	opts.seed = time(0);
	if (input.cmdOptionExists("-h") || argc < 2) {
		showHelpPage();
		return 0;
	}
	if (input.cmdOptionExists("-q")) opts.quadMetric = true;
	if (input.cmdOptionExists("-o1")) opts.adaptive = true;
	if (input.cmdOptionExists("-o2")) opts.anisotropic = true;
	if (input.cmdOptionExists("-o3")) {
		opts.anisotropic = true; 
		opts.quadMetric = true;
	}
	if (input.cmdOptionExists("-v")) opts.validation = true;
	if (input.cmdOptionExists("-n")) 
		opts.numCluster = input.getIntArg("-n", opts.numCluster);
	if (input.cmdOptionExists("-r")) 
		opts.ringLevel = input.getIntArg("-r", opts.ringLevel);
	if (input.cmdOptionExists("-seed")) 
		opts.seed = input.getIntArg("-seed", opts.seed);
	string in_str = input.getLastArg();
	string out_str = input.cmdOptionExists("-f") ? input.getCmdOption("-f") : genOutputFileName(in_str, opts);

	if (input.cmdOptionExists("-quad")) {
		MeshSimplify *ms = new MeshSimplify();
		FILE *_;
		_ = freopen(in_str.c_str(), "r", stdin);
		_ = freopen(out_str.c_str(), "w", stdout);
		ms->setRatio(input.getFloatArg("-ratio", 0.1));
		ms->input();
		ms->start();
		ms->output();
		delete ms;
		return 0;
	}
	VertexClustering* vc = new VertexClustering(opts);
	vc->input(in_str);
	vc->start();
	vc->build();
	vc->output(out_str);
	vc->stat();
	cout << "file saved in " << out_str << endl;
	return 0;
}

string VCFrontEnd::genOutputFileName(const string& in_file, const VertexClustering::Options& opts)
{
	stringstream ss;
	size_t k = in_file.find_last_of('/');
	string folder = ".", name = in_file;
	if (k != string::npos) {
		folder = in_file.substr(0, k);
		name = in_file.substr(k + 1);
	}
	name = name.substr(0, name.find_last_of('.'));
	string out_folder = folder + "/outputs";
#if defined(_WIN32)
	CreateDirectory(out_folder.c_str(), NULL);
#else
	system(("mkdir -p " + out_folder).c_str());
#endif
	ss << out_folder << "/" << name;
	if (opts.anisotropic)    ss << "_ais";
	else if (opts.adaptive)  ss << "_adp";
	else                     ss << "_uni";
	if (opts.quadMetric)     ss << "_q";
	if (opts.validation)     ss << "_v";
	if (opts.ringLevel != 2) ss << "_r" << opts.ringLevel;
	ss << "_" << opts.numCluster << ".obj";
	return ss.str();
}

void VCFrontEnd::showHelpPage()
{
	cout << "Usage: [options] <input_file>" << endl;
	cout << "-n INT    : set the number of clusters, default 400" << endl;
	cout << "-f PATH   : set the name of the output file, a default name will be generated if not provided" << endl;
	cout << "-r INT    : set the radius of the ring used by adaptive sampling, default 2" << endl;
	cout << "-seed INT : set the seed for srand(), default time(0)" << endl;
	cout << "-q        : use quadric based placement policy for vertex position" << endl;
	cout << "-o1       : enable adaptive sampling" << endl;
	cout << "-o2       : enable adaptive sampling with anisotropic metric" << endl;
	cout << "-o3       : shortcut for -o2 -q" << endl;
	cout << "-v        : enable verification for the clusters" << endl;
	cout << "-quad -f PATH -ratio FLOAT : run Quadric Error Metrics algorithm instead" << endl;
	cout << "-h        : show this help page" << endl;
}

