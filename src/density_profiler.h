#ifndef density_profiler_H
#define density_profiler_H
#include <string>
#include <vector>
#include <map>
using namespace std;
class gap_interval{
public:
	double start, stop;
	vector<double> X;
	vector<double> Y;
	gap_interval();
	gap_interval(double, double);
};

double get_table_mean_var(string, string,double, double, double);



#endif