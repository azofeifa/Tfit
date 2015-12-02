#ifndef template_H
#define template_H
#include "load.h"
#include <math.h>
#include <vector>
#include <iostream>
using namespace std;
vector<double> peak_bidirs(segment * );
int sample_centers(vector<double>, double);
void run_global_template_matching(vector<segment*>, 
	string, double, double, double, double, int, double, int );
void noise_global_template_matching(vector<segment*>, double);
extern double INF;
extern double  nINF;

#endif
