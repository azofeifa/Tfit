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
#include "model_selection.h"

using namespace std;
int main(int argc, char* argv[]){


	params * P 	= new params();
    P 			= readInParameters(argv);
    if (P->EXIT){
    	printf("exiting...\n");
    	return 0;
	}
    if (P->module=="MODEL"){


		//==========================================
		//Parameters
		int verbose 			= stoi(P->p["-v"]);
		string formatted_file 	= P->p["-i"]; 
		string out_file_dir 	= P->p["-o"] ;
		

		
		if (verbose){//show current user parameters...
			P->display();
		}
		//==========================================
		vector<segment*> segments	= load_EMGU_format_file(formatted_file, P->p["-chr"]);
		if (segments.empty()){
			cout<<"segments was not populated"<<endl;
			cout<<"exiting"<<endl;
			delete P;
			return 0;

		}
		BIN(segments, stod(P->p["-br"]), stod(P->p["-ns"]));
		//==========================================
		//Model Parameters
		
		clock_t t;
		chrono::time_point<chrono::system_clock> start, end;
		start = chrono::system_clock::now();

		t = clock();
		run_model_accross_segments(segments, P);
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
		cout<<"Loading Files"<<endl;
		map<string, vector<merged_interval*> > G 	= load_intervals(interval_file, pad); //load the forward and reverse strand intervals, merge accordingly...
		cout<<"Loaded Files"<<endl;
		//====================================================
		//making interval tree
		map<string, interval_tree *> A;
		typedef map<std::string, vector<merged_interval*>>::iterator it_type;
		cout<<"Making interval Tree"<<endl;
		for(it_type c = G.begin(); c != G.end(); c++) {
			A[c->first] 	= new interval_tree();
			A[c->first]->construct(c->second);
		}
		cout<<"Made interval Tree"<<endl;

		cout<<"Inserting Forward Strand"<<endl;
		insert_bedgraph(A, forward_strand_bedgraph_file, 1);
		cout<<"Inserted Forward Strand"<<endl;
		
		cout<<"Inserting Reverse Strand"<<endl;
		insert_bedgraph(A, reverse_strand_bedgraph_file, -1);
		cout<<"Inserted Reverse Strand"<<endl;
		
		cout<<"Writing OUT"<<endl;
		write_out(out_file_name, A);
		cout<<"DONE"<<endl;
	}
	else if(P->module=="SELECTION"){
		int verbose 			= stoi(P->p3["-v"]);
		
		if (verbose){//show current user parameters...
			P->display();
		}
		run_model_selection(P->p3["-i"], P->p3["-o"], stod(P->p3["-penality"]));
		
	}
		
	else {
		printf("Could not understand module or not provided...\n");
		printf("exiting\n");
		return 0;
	}
	delete P;
	return 1;
}