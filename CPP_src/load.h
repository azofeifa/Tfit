#ifndef load_H
#define load_H
#include <string>
#include <vector>
#include <map>

using namespace std;
class simple_c;
class model_component{
private:
	double mu, si, l, w_e, pi;
	double f_b, f_w;
	double r_a, r_w;
public:
	model_component(simple_c sc);
	model_component();
};

class all_model_components{
private:
	double ll;
	vector<model_component> all_components;
public:
	all_model_components();
	void insert_component(simple_c sc);
};



class bidir_preds{
private:
	double noise_ll ;
	map<int, all_model_components > G;
public:
	bidir_preds(double);
	bidir_preds();
	void insert_component(int, simple_c);
};




class classifier; //forward declare
class segment{
public:
	string chrom; 
	int start, stop;
	double minX, maxX;
	vector< vector<double> > forward;
	vector< vector<double> > reverse;
	int counts;
	segment(string, int , int);
	segment();
	string write_out();
	void bin(double, double, bool);
	void add(int, double, double);
	void add2(int, double, double);
	double N;
	double XN;
	double ** X;
	double SCALE;
	vector<vector<double> > bidirectional_bounds;
	vector<segment *> bidirectional_data;
	vector<int>  bidir_counts; //used for optimization of BIC?
	vector<int> bidirectional_N;
	void insert_bidirectional_data(int);

};




class interval{
public:
	string chrom;
	int start, stop, strand; //strand, 1 == forward, -1 == reverse
	interval();
	interval(string, int, int );
	vector<double> forward_x, forward_y, reverse_x, reverse_y;
	void insert(double, double, int);
	int hits;
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
	bool find(int, int);
	int get_hits(bool, int,int );
	void reset_hits();
	int get_total(int , int);


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
	bool find(int, int);
	int get_hits(bool, int, int);
	void reset_hits();
	int get_total(int, int);

};

vector<segment*> load_EMGU_format_file(string, string);

void BIN(vector<segment*>, int, double, bool);
map<string, vector<merged_interval*> > load_intervals(string, int);
void insert_bedgraph(map<string, interval_tree *>, string, int);
void write_out(string,map<string, interval_tree *> );

vector<segment*> load_bedgraphs_total(string, 
	string, int , double, string);
map<string, interval_tree *> load_bidir_bed_files(string,
	string);
#endif