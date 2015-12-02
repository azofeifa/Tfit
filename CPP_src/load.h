#ifndef load_H
#define load_H
#include <string>
#include <vector>
#include <map>
#include "read_in_parameters.h"
using namespace std;
class simple_c_free_mode;




class classifier; //forward declare



class segment{
public:
	string chrom; 
	int start, stop, ID, chrom_ID;
	double minX, maxX;
	vector< vector<double> > forward;
	vector< vector<double> > reverse;
	string strand ;
	int counts;
	vector<double> centers;
	vector<vector<double>> parameters; //for bootstraping
	map<int, vector<double> > variances;
	segment(string, int , int);
	segment(string, int , int, int);
	segment(string, int , int, int,string);
	segment();
	string write_out();
	void bin(double, double, bool);
	void add2(int, double, double);
	double N;
	double fN;
	double rN;
	double XN;
	double ** X;
	double SCALE;
	vector<vector<double> > bidirectional_bounds;
	vector<segment *> bidirectional_data;
	vector<int>  bidir_counts; //used for optimization of BIC?
	vector<int> bidirectional_N;
	vector<vector<double>> fitted_bidirs; //mu, si, l,pi
};

class node{
public:
	double center;
	int start, stop;
	node * left;
	node * right;
	vector<segment * > current;
	void retrieve_nodes(vector<segment * >&);
	void insert_coverage(vector<double>, int);
	node();
	node(vector<segment *>);
	void searchInterval(int, int, vector<int> &) ;
};

class segment_fits{
public:
	string chrom;
	int start, stop, TSS;
	double N;
	double N_pos, N_neg;
	map<int, double> M;
	map<int, string> parameters;
	int model;
	string ID;
	segment_fits();
	segment_fits(string, int, int, double, double, string);
	void get_model(double);
	string write();
};



namespace load{

	vector<segment_fits *> label_tss(string , vector<segment_fits *>   );

	void BIN(vector<segment*>, int, double, bool);


	vector<segment*> load_bedgraphs_total(string, 
		string, int , double, string,map<string, int>&,map<int, string>&);


	void write_out_bidirs(map<string , vector<vector<double> > >, string, string, int ,params *, int);

	vector<segment *> load_intervals_of_interest(string,map<int, string>&, params *);


	void collect_all_tmp_files(string , string, int, int );
	vector<segment* > insert_bedgraph_to_segment_joint(map<string, vector<segment *> >  , string , string , int);

	void write_out_models_from_free_mode(map<int, map<int, vector<simple_c_free_mode>  > >,
		params *,int,map<int, string>, int, string &);
	void clear_segments(vector<segment *> );

	vector<segment_fits *> load_K_models_out(string);
	void write_out_bidirectionals_ms_pen(vector<segment_fits*> , params * , int, int );

}

#endif