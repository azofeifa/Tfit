#ifndef across_segments_H
#define across_segments_H
#include "load.h"

void run_model_accross_segments(vector<segment*>, 
	int , int, int, double, 
	double, double, double, int, string, double);

void free_segments(vector<segment*>);


#endif