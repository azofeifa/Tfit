#ifndef across_segments_H
#define across_segments_H
#include "load.h"
#include "read_in_parameters.h"

string check_file(string, int);

void run_model_accross_segments(vector<segment*>, 
	params *);

void free_segments(vector<segment*>);


#endif