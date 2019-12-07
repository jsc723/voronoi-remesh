#ifndef CONFIG_H
#define CONFIG_H
// Stores some constants
class Config
{
public:
	Config(void);
	~Config(void);
	static const int MAX_NUM_VERTICES = 1000000; // maximum number of vertices
	static const int MAX_NUM_EDGES = 20000000; // maximum number of edges
	static const double EPS; // Negative infinity. Treated as zero
	static const double INF; // Positive infinity
};
#endif