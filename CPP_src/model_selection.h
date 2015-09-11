#ifndef model_selection_H
#define model_selection_H
#include <string>
#include <vector>
#include <map>
#include "load.h"
#include "model.h"
#include "MPI_comm.h"
using namespace std;
class model{
public:
	double ll;
	int K;
	vector<EMG> bidirs; 
	vector<UNI> elongations; 
	model(int, double);
	model();
	void add_component(bool, vector<string>, int);

};
class final_model_output{
public:
	string header;
	int K;
	double noise_ll;
	double k_ll;
	double scale;
	int start, stop, ID;
	string chrom;
	vector<rsimple_c> components;
	final_model_output(string, string, int,vector<rsimple_c>, double, double,double,int, int, int );
	final_model_output();
	string write_out_config();
	string write_out_bed();
	vector<vector<double> > get_bounds();

};


class S{
public:
	int start, stop;
	double N;
	string chrom;
	S(string , int , int, double);
	void add_model(int ,double);
	map<int, vector<model> > models;
	map<int, model> bests;
	model best;
	void pick_best(double );
	string print_best();
	string print_header();
	string getBestComponent(double );
	bool print_out;

};


void run_model_selection(string, string, double );
map<int, map<int, bidir_preds> > run_model_selection_bidir_template(
	map<int, map<int, bidir_preds> > , double );

vector<final_model_output> optimize_model_selection_bidirs(map<string, map<int, vector<rsimple_c> > >, params *);
vector<final_model_output> convert_to_final_model_output(map<string, map<int, vector<rsimple_c> > >, params *);
vector<final_model_output> convert_bidir_segs_to_final_model(map<string, vector<vector<double>>  > );


#endif