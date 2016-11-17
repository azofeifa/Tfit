#ifndef bootstrap_H
#define bootstrap_H
#include "load.h"
#include <vector>
#include "read_in_parameters.h"
#include <iostream>
#include <fstream>

using namespace std;


void run_bootstrap_across(vector<segment *>, params *, ofstream&  );

#endif