#include "load.h"
#include <iostream>
#include "across_segments.h"
#include <limits>
#include <math.h>
#include <errno.h>
using namespace std;
int main(){
	//==========================================
	//Parameters
	string formatted_file="/Users/joeyazo/Desktop/Lab/gro_seq_files/HCT116/EMG_out_files/test_file_2.tsv"; 
	string spec_all 	= "chrN";
	int bin_resolution 	= 300;
	double scale 		= 100;
	//==========================================
	vector<segment*> segments	= load_EMGU_format_file(formatted_file, spec_all);
	if (segments.empty()){
		cout<<"segments was not populated"<<endl;
		cout<<"exiting"<<endl;

	}
	BIN(segments, bin_resolution, scale);
	//==========================================
	//Model Parameters
	int maxK 	= 3; //max number of models to try
	int rounds 	= 10; //number of random seeds


	run_model_accross_segments(segments, maxK, rounds);
	
	
	return 1;
}