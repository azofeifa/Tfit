#include "load.h"
#include <iostream>
#include "across_segments.h"
#include <limits>
#include <math.h>
#include <errno.h>
#include <time.h>  
#include <stdio.h>   
#include <chrono>
#include "read_in_parameters.h"
using namespace std;
int main(int argc, char* argv[]){


	params * P 	= new params();
    P 			= readInParameters(argv);




	//==========================================
	//Parameters
	int verbose 			= stoi(P->p["-v"]);
	string formatted_file 	= P->p["-i"]; 
	string out_file_dir 	= P->p["-o"] ;
	string spec_all 		= P->p["-chr"];
	int bin_resolution 		= stoi(P->p["-br"]);
	double scale 			= stod(P->p["-ns"]);
	int num_proc 			= stoi(P->p["-np"]);
	double move 			= stod(P->p["-move"]);
	double max_noise 		= stod(P->p["-max_noise"]);
	double convergence_tresh= stod(P->p["-ct"]);
	int max_iterations 		= stoi(P->p["-mi"]);
	bool pp_r 				= bool(stoi(P->p["-pp_r"]));
	bool pp_k 				= bool(stoi(P->p["-pp_k"]));
	double r_mu 			= stod(P->p["-r_mu"]);
	
	if (verbose){//show current user parameters...
		P->display();
	}
	//==========================================
	vector<segment*> segments	= load_EMGU_format_file(formatted_file, spec_all);
	if (segments.empty()){
		cout<<"segments was not populated"<<endl;
		cout<<"exiting"<<endl;
		delete P;
		return 0;

	}
	BIN(segments, bin_resolution, scale);
	//==========================================
	//Model Parameters
	int maxK 	= stoi(P->p["-maxK"]); //max number of models to try
	int rounds 	= stoi(P->p["-rounds"]); //number of random seeds

	clock_t t;
	chrono::time_point<chrono::system_clock> start, end;
	start = chrono::system_clock::now();

	t = clock();
	run_model_accross_segments(segments, maxK, rounds, num_proc, scale, move, 
		max_noise, convergence_tresh, max_iterations,out_file_dir, pp_r, pp_k,
		r_mu);
	free_segments(segments);
	end = chrono::system_clock::now();
	int elapsed_seconds = std::chrono::duration_cast<std::chrono::milliseconds >
                             (end-start).count();
	t = clock() - t;
	printf ("CPU: %d , %f\n",t,((float)t)/CLOCKS_PER_SEC);
	printf("Wall Time: %f\n", (elapsed_seconds/1000.));
	delete P;
	return 1;
}