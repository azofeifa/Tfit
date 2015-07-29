#include "load.h"
#include "across_segments.h"
#include "model.h"
#include "template_matching.h"
#include <iostream>
#include "omp.h"
#include <fstream>
#include <map>
#include <time.h>
using namespace std;

const string currentDateTime() {
	time_t     now = time(0);
	struct tm  tstruct;
	char       buf[80];
	tstruct = *localtime(&now);
	strftime(buf, sizeof(buf), "%Y-%m-%d,%X", &tstruct);

	return buf;
}

string check_file(string FILE, int i){ //don't want to write over an existing file
	string template_file 	= FILE + to_string(i);
	ifstream FH(template_file);
	if (FH){
		FH.close();
		return	check_file(FILE, i+1);
	}
	return template_file;
}

map<int,vector<classifier> > initialize_data_struct(int minK, int maxK, 
	int rounds, int num_proc, double scale, 
	double max_noise,  double convergence_tresh, 
	int max_iterations, double r_mu, double ALPHA_0, double BETA_0, 
	double ALPHA_1, double BETA_1, double ALPHA_2, double ALPHA_3){
	map<int,vector<classifier> >  data_struct;	
	for (int k = minK; k<=maxK;k++){
		data_struct[k] 	= vector<classifier>(rounds);
		for (int r = 0; r<rounds; r++){
			data_struct[k][r] 	= classifier(k,  convergence_tresh, max_iterations, 
					max_noise,r_mu, ALPHA_0, BETA_0, ALPHA_1, BETA_1, ALPHA_2, ALPHA_3);
		}
	}
	return data_struct;
}

map<int, classifier> getMax(map<int,vector<classifier> > DS){
	typedef map<int, vector<classifier>>::iterator it_type;
	map<int, classifier> reduced;
	double MAX=nINF;
	for(it_type K = DS.begin(); K != DS.end(); K++) {
		for (int i = 0; i < K->second.size(); i ++){
			if (DS[K->first][i].ll > MAX and DS[K->first][i].converged) {
				reduced[K->first] 	= DS[K->first][i];
				MAX 				= DS[K->first][i].ll;
			}
		}
	}
	return reduced;

}



void run_model_accross_segments(vector<segment*> segments, params *P){
	typedef map<int, classifier>::iterator it;
	



	int N 			= segments.size();
	string out_dir 	= P->p["-o"]; 
	string out_file_template 	= out_dir+"model_fits_out_" + P->p["-chr"]+"_";
	string out_file 			= check_file(out_file_template, 1);
	
	

	ofstream FHW;
	
	FHW.open(out_file);

	string header 	= currentDateTime() + ",";
	header+="minK: " + P->p["-minK"] + ",";
	header+="maxK: " +  P->p["-maxK"] + ",";
	header+="scale: " +  P->p["-ns"] + ",";
	header+="bin_resolution: " + P->p["-br"] + ",";
	header+="convergence_tresh: " + P->p["-ct"] + ",";
	header+="max_iterations: " + P->p["-mi"] + ",";
	header+="spec_chrom: " + P->p["-chr"] + ",";
	header+="max_noise: " + P->p["-max_noise"] + ",";
	header+="r_mu: " + P->p["-r_mu"] +",";
	header+="rounds: " + P->p["-rounds"]+",";
	header+="ALPHA_0: " + P->p["-ALPHA_0"] +",";
	header+="BETA_0: " + P->p["-BETA_0"] +",";
	header+="ALPHA_1: " + P->p["-ALPHA_1"] +",";
	header+="BETA_1: " + P->p["-BETA_1"] +",";
	header+="ALPHA_2: " + P->p["-ALPHA_2"] +",";
	header+="ALPHA_3: " + P->p["-ALPHA_3"]  ;
	
	
	
	
	header+="\n";
	FHW<<header;

	int minK 					= stoi(P->p["-minK"]);
	int maxK 					= stoi(P->p["-maxK"]);
	double max_noise 			= stod(P->p["-max_noise"]);
	double convergence_tresh 	= stod(P->p["-ct"]);
	int max_iterations 			= stoi(P->p["-mi"]);
	double r_mu 				= stod(P->p["-r_mu"]);
	int num_proc 				= stoi(P->p["-np"]);
	int rounds 					= stoi(P->p["-rounds"]);
	double scale 				= stod(P->p["-ns"]);
	double ALPHA_0 				= stod(P->p["-ALPHA_0"]);
	double BETA_0 				= stod(P->p["-BETA_0"]);
	double ALPHA_1 				= stod(P->p["-ALPHA_1"]);
	double BETA_1 				= stod(P->p["-BETA_1"]);
	double ALPHA_2 				= stod(P->p["-ALPHA_2"]);
	double ALPHA_3 				= stod(P->p["-ALPHA_3"]);
	
	for (int i = 0; i < N; i++ ){
		if (segments[i]->N > 0){
			vector<double> mu_seeds 		=  peak_bidirs(segments[i]);
			map<int,vector<classifier> > DS = initialize_data_struct(minK, maxK, 
				rounds, num_proc, scale, max_noise,  
				convergence_tresh, max_iterations,r_mu, ALPHA_0,  BETA_0, 
	 		ALPHA_1,  BETA_1,  ALPHA_2,  ALPHA_3);
			

			classifier clf(0, convergence_tresh, max_iterations, 
						max_noise,r_mu, ALPHA_0, BETA_0, ALPHA_1, BETA_1, ALPHA_2, ALPHA_3);
			clf.fit(segments[i], mu_seeds);
			
			FHW<<segments[i]->write_out();
			FHW<<"~0"<<","<<to_string(clf.ll)<<",1,0"<<endl;
			FHW<<"U: "<<to_string(segments[i]->minX)<<","<<to_string(segments[i]->maxX)<<",1,"<<to_string(clf.pi)<<endl;
			for (int k = minK; k <=maxK;k++ ){
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
			for (int k = minK; k<=maxK; k++){
				if (reduced.find(k) !=reduced.end()){
					classifier K 	= reduced[k];
					FHW<<"~"<<to_string(k)<<","<<to_string(K.ll)<<","<<to_string(K.converged)<<","<<to_string(K.last_diff)<<endl;
					FHW<<K.print_out_components();
				}else{
					FHW<<"~"<<to_string(k)<<","<<to_string(nINF)<<","<<to_string(0)<<","<<to_string(INF)<<endl;
				}
			}
			
		}		
	}
}

void free_segments(vector<segment*> segments){
	for (int i = 0; i < segments.size(); i++){
		delete segments[i];
	}
}


