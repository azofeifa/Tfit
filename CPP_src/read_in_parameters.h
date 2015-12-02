#ifndef read_in_parameters_H
#define read_in_parameters_H

#include <iostream>
#include <fstream>
#include <unistd.h>
#include <string>
#include <vector>
#include <map>
using namespace std;


class params{
public:
	map<string, string> p;
	params();
	void display(int,int);
	void help();
	int N;
	string module;
	bool bidir;
	bool model;
	bool select;

	map<string, string> p2;

	map<string, string> p3;
	map<string, string> p4;
	map<string, string> p5;
	map<string, string> p6;
	
	char * isIntGroup[8] = {"-pad", "-minK", "-maxK", 
						 "-rounds", "-mi", "-MLE", "-elon", "-merge"};

	char * isDecGroup[12]  = {  "-br","-ns", "-ct",
						"-max_noise",    "-r_mu",
						"-ALPHA_0", "-ALPHA_1", "-ALPHA_2", "-BETA_0", "-BETA_1",
						"-bct", "-ms_pen"   };  
	char * isPathGroup[8] = {"-config", "-i", "-j", "-k", "-tss", "-log_out", "-o", "-q"};
	string get_header(int);
	
	vector<string> validate_parameters();
	bool EXIT;
	
};

void fillInOptions(char*,params);

params * readInParameters(char**);

const std::string currentDateTime();
void fill_in_bidir_boostrap(params *);
int read_in_parameters( char**, params *, int );
#endif