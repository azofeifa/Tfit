#include "load.h"
#include "across_segments.h"
#include "model.h"
#include "template_matching.h"
#include <iostream>
#include "omp.h"
#include <fstream>
#include <map>
using namespace std;
string check_file(string FILE, int i){ //don't want to write over an existing file
	string template_file 	= FILE + to_string(i);
	ifstream FH(template_file);
	if (FH){
		FH.close();
		return	check_file(FILE, i+1);
	}
	return template_file;
}

map<int,vector<classifier> > initialize_data_struct(int maxK, 
	int rounds, int num_proc, double scale, double move, 
	double max_noise,  double convergence_tresh, 
	int max_iterations, double r_mu){
	map<int,vector<classifier> >  data_struct;	
	for (int k = 1; k<=maxK;k++){
		data_struct[k] 	= vector<classifier>(rounds);
		for (int r = 0; r<rounds; r++){
			data_struct[k][r] 	= classifier(k,  convergence_tresh, max_iterations, 
					max_noise, move,r_mu);
		}
	}
	return data_struct;
}

map<int, classifier> getMax(map<int,vector<classifier> > DS){
	typedef map<int, vector<classifier>>::iterator it_type;
	map<int, classifier> reduced;
	double MAX=0;
	for(it_type K = DS.begin(); K != DS.end(); K++) {
		MAX=0;
		for (int i = 0; i < K->second.size(); i ++){
			if (MAX==0){
				reduced[K->first] 	= DS[K->first][i];
				MAX 				= DS[K->first][i].ll;
			}else if (DS[K->first][i].ll > MAX) {
				reduced[K->first] 	= DS[K->first][i];
				MAX 				= DS[K->first][i].ll;
			}
		}
	}
	return reduced;

}



void run_model_accross_segments(vector<segment*> segments, 
	int maxK, int rounds, int num_proc, double scale, double move, 
		double max_noise,  double convergence_tresh, int max_iterations,
		string out_dir, double r_mu){
	typedef map<int, classifier>::iterator it;
	
	int N 	= segments.size();
	string out_file_template 	= out_dir+"model_fits_out_";
	string out_file 			= check_file(out_file_template, 1);
	
	


	ofstream FHW;
	
	FHW.open(out_file);
	for (int i = 0; i < N; i++ ){
	 	vector<double> mu_seeds 		=  peak_bidirs(segments[i]);
	 			
		map<int,vector<classifier> > DS = initialize_data_struct(maxK, 
			rounds, num_proc, scale,  move, max_noise,  
			convergence_tresh, max_iterations,r_mu);
		FHW<<segments[i]->write_out();
		classifier clf(0, convergence_tresh, max_iterations, 
					max_noise, move,r_mu);
		clf.fit(segments[i], mu_seeds);

		FHW<<"~0"<<","<<to_string(clf.ll)<<",1,0"<<endl;
		FHW<<"U: "<<to_string(segments[i]->minX)<<","<<to_string(segments[i]->maxX)<<",1,"<<to_string(clf.pi)<<endl;

		//#pragma omp parallel for num_threads(num_proc_k)
		//#pragma omp parallel for num_threads(num_proc) collapse(2)
		for (int k = 1; k <=maxK;k++ ){
			#pragma omp parallel for num_threads(num_proc)
			for (int j = 0; j < rounds; j++){
				DS[k][j].fit(segments[i], mu_seeds);
			}
		}
		if (segments[i]->N > 0){
			vector<double> mu_seeds 		=  peak_bidirs(segments[i]);
			map<int,vector<classifier> > DS = initialize_data_struct(maxK, 
				rounds, num_proc, scale,  move, max_noise,  
				convergence_tresh, max_iterations,r_mu);
			FHW<<segments[i]->write_out();
			classifier clf(0, convergence_tresh, max_iterations, 
						max_noise, move,r_mu);
			clf.fit(segments[i], mu_seeds);
			
			FHW<<"~0"<<","<<to_string(clf.ll)<<",1,0"<<endl;
			FHW<<"U: "<<to_string(segments[i]->minX)<<","<<to_string(segments[i]->maxX)<<",1,"<<to_string(clf.pi)<<endl;

			for (int k = 1; k <=maxK;k++ ){
				#pragma omp parallel for num_threads(num_proc)
				for (int j = 0; j < rounds; j++){
					DS[k][j].fit(segments[i], mu_seeds);
				}
			}
			//*********************************************
			//we want to output for each model, k, the highest log likelihood estimate
			//and then the respective parameter estimates for each component
			map<int, classifier> reduced 	= getMax(DS);
			//okay now lets write it out to the file
			for(it  K = reduced.begin(); K != reduced.end(); K++) {
				FHW<<"~"<<to_string(K->first)<<","<<to_string(K->second.ll)<<","<<to_string(K->second.converged)<<","<<to_string(K->second.last_diff)<<endl;
				FHW<<K->second.print_out_components();
			}
		}		
	}
}

void free_segments(vector<segment*> segments){
	for (int i = 0; i < segments.size(); i++){
		delete segments[i];
	}
}


