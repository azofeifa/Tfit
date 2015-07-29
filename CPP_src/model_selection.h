#ifndef model_selection_H
#define model_selection_H
#include <string>
#include <vector>

#include "model.h"
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

#endif