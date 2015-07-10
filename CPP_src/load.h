#ifndef load_H
#define load_H

#include <string>
#include <vector>
using namespace std;
class segment{
public:
	string chrom; 
	int start, stop;
	double minX, maxX;
	vector< vector<double> > forward;
	vector< vector<double> > reverse;
	segment(string, int , int);
	segment();
	void bin(int, double);
	void add(int, double, double);
	double N;
	double XN;
	double ** X;
	double SCALE;
};

vector<segment*> load_EMGU_format_file(string, string);
void BIN(vector<segment*>, int, double);
#endif