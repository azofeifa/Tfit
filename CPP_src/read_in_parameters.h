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
	void display();
	void help();
	int N;
	string module;
	map<string, string> p2;

	map<string, string> p3;
	bool EXIT;
	
};

void fillInOptions(char*,params);

params * readInParameters(char**);

#endif