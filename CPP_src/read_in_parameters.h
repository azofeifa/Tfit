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
	void display(int);
	void help();
	int N;
	string module;
	map<string, string> p2;

	map<string, string> p3;
	map<string, string> p4;
	map<string, string> p5;
	string get_header(int);
	
	bool EXIT;
	
};

void fillInOptions(char*,params);

params * readInParameters(char**);

const std::string currentDateTime();
#endif