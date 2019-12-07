#ifndef VCFRONTEND_H
#define VCFRONTEND_H
#include "vertexClustering.h"
using std::vector;
using std::string;

namespace VCTools {
	//ref: https://stackoverflow.com/questions/865668/how-to-parse-command-line-arguments-in-c
	class InputParser {
	public:
		InputParser(int argc, char** argv) {
			for (int i = 1; i < argc; ++i)
				this->tokens.push_back(string(argv[i]));
		}
		/// @author iain
		string getCmdOption(const string & option) const {
			vector<string>::const_iterator itr;
			itr = find(this->tokens.begin(), this->tokens.end(), option);
			if (itr != this->tokens.end() && ++itr != this->tokens.end()) {
				return *itr;
			}
			static const string empty_string("");
			return empty_string;
		}
		/// @author iain
		bool cmdOptionExists(const string & option) const {
			return find(this->tokens.begin(), this->tokens.end(), option)
				!= this->tokens.end();
		}
		string getLastArg() const {
			return *tokens.rbegin();
		}
		int getIntArg(const string& arg, int def = -1) {
			string nr = getCmdOption(arg);
			if (nr.size() > 0 && std::all_of(nr.begin(), nr.end(), ::isdigit))
				return atoi(nr.c_str());
			return def;
		}
		double getFloatArg(const string& arg, double def = -1) {
			string nr = getCmdOption(arg);
			if (nr.size() > 0 && std::all_of(nr.begin(), nr.end(), 
				[](unsigned char c) {return ::isdigit(c) || c == '.'; }))
				return atof(nr.c_str());
			return def;
		}
	private:
		vector <string> tokens;
	};
}


class VCFrontEnd {
public:
	static int exec(int argc, char** argv);
	static string genOutputFileName(const string &in_file, const VertexClustering::Options& opts);
	static void showHelpPage();
};

#endif // !VCFRONTEND_H
