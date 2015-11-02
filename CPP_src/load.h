#ifndef load_H
#define load_H
#include <string>
#include <vector>
#include <map>
#include "read_in_parameters.h"
using namespace std;
class simple_c;
class simple_c_free_mode;
class final_model_output;
struct single_simple_c;
struct boostrap_struct;
class model_component{
public:
	double mu, si, l, w_e, pi;
	double f_b, f_w;
	double r_a, r_w;

	model_component(simple_c sc);
	model_component();
};

class all_model_components{
public:
	double ll;
	vector<model_component> all_components;
	all_model_components();
	void insert_component(simple_c sc, double);
};



class bidir_preds{
public:
	double noise_ll ;
	double N;
	map<int, all_model_components > G;
	all_model_components best_component;
	bidir_preds(double, double);
	bidir_preds();
	void insert_component(int, simple_c,double);
	void model_selection(double);
	int argK;
	double BIC_score;
	all_model_components arg_all_model_components;
};




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
	vector<vector<double>> fitted_bidirs; //mu, si, l,pi
	void insert_bidirectional_data(int);
	void add_fitted_bidir(vector<double>);
};




class interval{
public:
	string chrom;
	int ID;
	int counts;
	int start, stop, strand; //strand, 1 == forward, -1 == reverse
	bool EMPTY;
	string STRAND;
	vector<vector<double>> parameters; //for bootstraping

	interval();
	interval(string, int, int );
	interval(string, int, int , int );
	interval(string, int, int, int , string,vector<vector<double>>);
	interval(string, int, int, int , string,vector<vector<double>>, int);
	
	vector<double> forward_x, forward_y, reverse_x, reverse_y;
	void insert(double, double, int);
	int hits;
};

class merged_interval{
public:
	string chrom;
	int start, stop;
	int id;
	vector<vector<double>> parameters;
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
	interval get_interval(int, int);


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
	interval get_interval(int, int);
};

vector<segment*> load_EMGU_format_file(string, string);

void BIN(vector<segment*>, int, double, bool);
map<string, vector<merged_interval*> > load_intervals(string, int);
void insert_bedgraph(map<string, interval_tree *>, string, int);
void write_out(string,map<string, interval_tree *> );

vector<segment*> load_bedgraphs_total(string, 
	string, int , double, string,map<string, int>&,map<int, string>&);
vector<segment*> load_bedgraphs_single(string, int , double, string);
map<string, interval_tree *> load_bidir_bed_files(string,
	string);

void write_out_bidir_fits( vector<segment*>, 
	map<int, map<int, bidir_preds> >, params *);

void write_out_bidirs(map<string , vector<vector<double> > >, string, string, int ,params *);

vector<segment *> bidir_to_segment(map<string , vector<vector<double> > >, 
	string , string, int, string, map<string, int> );
vector<segment *>  combind_bidir_fits_with_intervals_of_interest(vector<final_model_output> , vector<segment *> );

void write_out_MLE_model_info(vector<final_model_output>, params *, string, int);

vector<segment *> load_intervals_of_interest(string,map<int, string>&, int, string);

vector<segment *> insert_bedgraph_to_segment(map<string, vector<segment *> > , string, string, int);

void write_gtf_file_model_fits(vector<final_model_output>, params *);

vector<segment* > insert_bedgraph_to_segment_single(map<string, vector<segment *> > , string, int);

void write_out_single_simple_c(vector<single_simple_c>, map<int, string> , params * );
void write_out_single_simple_c_marks(vector<single_simple_c>, map<int, string> , params * );

void write_config_file_model_fits(vector<final_model_output> , map<int, string>, params * );


void collect_all_tmp_files(string , string, int, int );
vector<segment* > insert_bedgraph_to_segment_joint(map<string, vector<segment *> >  , string , string , int);

void get_noise_mean_var(string, string, double *, double *);

void write_out_models_from_free_mode(map<int, map<int, vector<simple_c_free_mode>  > >,params *,int,map<int, string>);

map<string, vector<segment *> > load_bidir_predictions( params *,
	vector<int>, map<string, int>&, map<int, string>&  );
vector<vector<int> > get_line_start_stops(params * , int );
void write_bootstrap(vector<boostrap_struct> , map<int, string> , params * ,int );
vector<segment *> merge_intervals_of_interest(vector<segment *>);

#endif