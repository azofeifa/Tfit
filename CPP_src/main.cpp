#include "load.h"
#include <iostream>
#include "across_segments.h"
#include <limits>
#include <math.h>
#include <errno.h>
#include <time.h>  
#include <stdio.h>   
#include <chrono>
#include <map>
#include "read_in_parameters.h"

using namespace std;
int main(int argc, char* argv[]){


	params * P 	= new params();
    P 			= readInParameters(argv);

    if (P->module=="MODEL"){


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
		int minK 	= stoi(P->p["-minK"]);
		int rounds 	= stoi(P->p["-rounds"]); //number of random seeds

		clock_t t;
		chrono::time_point<chrono::system_clock> start, end;
		start = chrono::system_clock::now();

		t = clock();
		run_model_accross_segments(segments, minK, maxK, rounds, num_proc, scale, move, 
			max_noise, convergence_tresh, max_iterations,out_file_dir, 
			r_mu,spec_all);
		free_segments(segments);
		end = chrono::system_clock::now();
		int elapsed_seconds = std::chrono::duration_cast<std::chrono::milliseconds >
	                             (end-start).count();
		t = clock() - t;
		printf ("CPU: %d , %f\n",t,((float)t)/CLOCKS_PER_SEC);
		printf("Wall Time: %f\n", (elapsed_seconds/1000.));
	}else if (P->module=="FORMAT"){
		//==================================================================
		//Necessary FILES
		string interval_file 					= P->p2["-i"];
		string forward_strand_bedgraph_file 	= P->p2["-j"];
		string reverse_strand_bedgraph_file 	= P->p2["-k"];
		
		string out_file_name 					= P->p2["-o"];
		//==================================================================
		//Some cosmetic parameters
		int pad 								= stoi(P->p2["-pad"]);
		bool verbose 							= bool(stoi(P->p2["-v"]));
		//==================================================================
		if (verbose){//show current user parameters...
			P->display();
		}
		map<string, vector<merged_interval*> > G 	= load_intervals(interval_file); //load the forward and reverse strand intervals, merge accordingly...
		//====================================================
		//making interval tree
		map<string, interval_tree *> A;
		typedef map<std::string, vector<merged_interval*>>::iterator it_type;
		for(it_type c = G.begin(); c != G.end(); c++) {
			A[c->first] 	= new interval_tree();
			A[c->first]->construct(c->second);
		}
		insert_bedgraph(A, forward_strand_bedgraph_file, 1);
		
		insert_bedgraph(A, reverse_strand_bedgraph_file, -1);

		write_out(out_file_name, A);

	





	}else {
		printf("Could not understand module or not provided...\n");
		printf("exiting\n");
		return 0;
	}
	delete P;
	return 1;
}