#include "load.h"
#include "across_segments.h"
#include "model.h"
#include "template_matching.h"
#include <iostream>
using namespace std;
void run_model_accross_segments(vector<segment*> segments, int maxK, int rounds){
	int N 	= segments.size();
	for (int i = 0; i < N; i++ ){
		vector<double> mu_seeds 	=  peak_bidirs(segments[i]);
		
		for (int k = 1; k <=maxK;k++ ){
			for (int j = 0; j < rounds; j++){
				classifier clf(k);
				clf.fit(segments[i], mu_seeds);
			}
		}
			
	}
}