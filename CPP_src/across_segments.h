#ifndef across_segments_H
#define across_segments_H
#include "load.h"
#include "read_in_parameters.h"
#include "model.h"

string check_file(string, int);

void run_model_accross_segments(vector<segment*>, 
	params *);
void free_segments(vector<segment*>);

struct simple_c{
	double ll; //loglikelihood of bidirectional(s) 
	
	double noise_ll; //loglikelihood of noise component
	int IDS[4];

	//IDS[0] 	= segment that it belongs to...
	//IDS[1] 	= complexity, number of fitted models
	//IDS[2] 	= number of predicted bidirs
	//IDS[3] 	= bidir prediction, (possible merged segment)
	double ps[9];
	
//	string write_out();

};
vector<simple_c> run_model_accross_segments_template(vector<segment*>, 
	params *);

#endif