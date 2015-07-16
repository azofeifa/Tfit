#ifndef load_H
#define load_H

#include <string>
#include <vector>
#include <map>
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
	string write_out();
	void bin(int, double);
	void add(int, double, double);
	double N;
	double XN;
	double ** X;
	double SCALE;

};




class interval{
public:
	string chrom;
	int start, stop, strand; //strand, 1 == forward, -1 == reverse
	interval();
	interval(string, int, int );
	vector<double> forward_x, forward_y, reverse_x, reverse_y;
	void insert(double, double, int);
};

class merged_interval{
public:
	string chrom;
	int start, stop;
	int id;
	merged_interval * left;
	merged_interval * right;
	vector<interval> intervals; 
	merged_interval();
	merged_interval(int, int, interval, int);
	void update(interval);
	void insert(double , double , int );


};

class interval_tree{
public:
	interval_tree();
	merged_interval * current;
	interval_tree * left;
	interval_tree * right;
	int ID;
	void insert(double , double , int);
	void construct(vector<merged_interval*>);
	int getDepth();
	void insert_into_array(merged_interval *, int);
};

vector<segment*> load_EMGU_format_file(string, string);

void BIN(vector<segment*>, int, double);
map<string, vector<merged_interval*> > load_intervals(string);
void insert_bedgraph(map<string, interval_tree *>, string, int);
void write_out(string,map<string, interval_tree *> );
#endif