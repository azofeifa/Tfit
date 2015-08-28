#include <mpi.h>
#include "load.h"
#include "model.h"
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
#include <thread>
#include "template_matching.h"
#include "MPI_comm.h"
using namespace std;


int main(int argc, char* argv[]){
	MPI::Init(argc, argv);
	params * P 	= new params();
    P 			= readInParameters(argv);
    if (P->EXIT){
    	printf("exiting...\n");
    	MPI::Finalize();
    	return 1;
	}
	if (P->module == "BIDIR"){
		int nprocs = MPI::COMM_WORLD.Get_size();
		int rank = MPI::COMM_WORLD.Get_rank();
	    
		int verbose 			= stoi(P->p4["-v"]);
		if (verbose and rank==0){//show current user parameters...
			P->display();
			printf("Running on %d node(s)\n", nprocs );
		}

		//load bed graph files into map<string, double **>;
		string forward_bedgraph 	= P->p4["-i"] ;
		string reverse_bedgraph 	= P->p4["-j"] ;
		string optimize_directory 	= P->p4["-k"];
		if (P->p4["-optimize"]=="0"){
			optimize_directory 		= "";
		}

		string out_file_dir 		= P->p4["-o"] ;
		int BINS 					= stoi(P->p4["-br"]);
		double scale 				= stod(P->p4["-ns"]);
		double window 				= stod(P->p4["-window"])/ scale;
		double ct 					= stod(P->p4["-ct"]);
		double density 				= stod(P->p4["-density"]) / scale;
		int opt_res 				= stod(P->p4["-opt_res"]);
		int np 						= stoi(P->p4["-np"]);
		string spec_chrom 			= P->p4["-chr"];
		
		vector<segment*> segments 	= load_bedgraphs_total(forward_bedgraph, 
			reverse_bedgraph, BINS, scale, spec_chrom);
		if (segments.empty()){
			printf("segments not populated, exiting...\n");
			MPI::Finalize();
			return 1;
		}
		vector<segment*> all_segments  	= segments;
		
		segments 						= slice_segments(segments, rank, nprocs);
		if (not optimize_directory.empty() and rank==0 and not segments.empty() ){
			map<string, interval_tree *> I = load_bidir_bed_files(optimize_directory, 
				spec_chrom);
			if (I.empty()){
				printf("exiting...\n");
				MPI::Finalize();
				return 1;
			}
			clock_t t;
			chrono::time_point<chrono::system_clock> start, end;
			start = chrono::system_clock::now();

			t = clock();
			optimize(I, segments, scale, opt_res, out_file_dir, spec_chrom, np);
			end = chrono::system_clock::now();
			int elapsed_seconds = std::chrono::duration_cast<std::chrono::milliseconds >
		                             (end-start).count();
			t = clock() - t;
			printf("CPU: %f, %f\n",t,((float)t)/CLOCKS_PER_SEC);
			printf("Wall Time: %f\n", (elapsed_seconds/1000.));
			free_segments(segments);
			MPI::Finalize();
			return 0;			
		}else if(optimize_directory.empty() and not segments.empty() ){
			run_global_template_matching(segments, out_file_dir, window, 
				density,scale,ct, np,-0.1 );	
		}
		map<string , vector<vector<double> > > G = gather_all_bidir_predicitions(all_segments, segments , rank, nprocs, out_file_dir);
		vector<segment *> bidir_segments;
		if (not G.empty()  ){
			bidir_segments 	= bidir_to_segment( G, forward_bedgraph,reverse_bedgraph, stoi(P->p4["-pad"]));
		}
		vector<simple_c> fits;
		if (not bidir_segments.empty()){
			BIN(bidir_segments, stod(P->p4["-br"]), stod(P->p4["-ns"]),true );
			fits 			= run_model_accross_segments_to_simple_c(bidir_segments, P);
		}
		map<string, map<int, vector<rsimple_c> > > rcG 	= gather_all_simple_c_fits(bidir_segments, fits, rank, nprocs);
		if (rank==0 and not rcG.empty() ){//perform and optimize model selection based on number of bidir counts
			vector<final_model_output> 	A  				= optimize_model_selection_bidirs(rcG, P);
		}


		if (not segments.empty()){
			free_segments(segments);
		}
		MPI::Finalize();
	
		return 0;
	}


    else if (P->module=="MODEL"){
    	

    	//==========================================
		//Parameters
		int verbose 			= stoi(P->p["-v"]);
		int template_match 		= stoi(P->p["-template"]);
		string formatted_file 	= P->p["-i"]; 
		string out_file_dir 	= P->p["-o"] ;
		
	    // get the number of processes, and the id of this process
	    int rank = MPI::COMM_WORLD.Get_rank();
	    int nprocs = MPI::COMM_WORLD.Get_size();

		bool root 	= (rank==0)	;
		
		if (verbose and rank==0){//show current user parameters...
			P->display();
			printf("Running on %d node(s)\n", nprocs );
		}

		//==========================================
		vector<segment*> segments		= load_EMGU_format_file(formatted_file, P->p["-chr"]);
		vector<segment*> all_segments  	= segments;
		if (segments.empty()){
			cout<<"segments was not populated"<<endl;
			cout<<"exiting"<<endl;
			delete P;
			MPI::Finalize();
			return 1;

		}
		segments 	= slice_segments(segments, rank, nprocs);
		int interval_number;
		map<int, int> bidir_number_table;
		
		if (not template_match){
			BIN(segments, stod(P->p["-br"]), stod(P->p["-ns"]), true);
		}else{
			BIN(segments, stod(P->p["-br"]), stod(P->p["-ns"]), false );
			//need to write a function to insert eRNA predictions 
			//for each segment
			run_global_template_matching(segments, "", 
				stod(P->p["-window"])/stod(P->p["-ns"]), stod(P->p["-density"])/stod(P->p["-ns"]) ,
				stod(P->p["-ns"]) ,stod(P->p["-ct"]), 
				stoi(P->p["-np"]) ,0);
			//now we want to collapse all segments into one large vector<segment *>
		}
		
		//==========================================
		//Model Parameters
		clock_t t;
		chrono::time_point<chrono::system_clock> start, end;
		start = chrono::system_clock::now();

		t = clock();
		vector<simple_c> fits;
		vector<simple_c> all_fits;
		map<int, map<int, bidir_preds> > G;
		if (not template_match){
			run_model_accross_segments(segments, P);
		}else{//need to write a function to fit bidirectionals in isolation...
			fits 	= run_model_accross_segments_template(segments, P);
		}
		if (root){
			bidir_number_table 	= get_all_bidir_sizes(fits, nprocs);
			typedef map<int, int>::iterator it_type;
			
			for (it_type cc = bidir_number_table.begin(); cc != bidir_number_table.end() ; cc++){
				printf("Node %d ran %d model fits\n", cc->first, cc->second );
			}
		}else{
			send_bidir_size(fits);
		}
		//now we want to collate all the results;
		if (root){
			G 	= gather_all_simple_c_fits(fits, bidir_number_table, interval_number, nprocs);
		}else{
			send_all_simple_c_fits(fits);
		}
		//need to perform model selection
		if (root){
			G 	= run_model_selection_bidir_template(G, 1.);
			write_out_bidir_fits(all_segments, G, P);
		}
		free_segments(segments);
		end = chrono::system_clock::now();
		int elapsed_seconds = std::chrono::duration_cast<std::chrono::milliseconds >
	                             (end-start).count();
		t = clock() - t;
		if (root){
			printf("Wall Time %f\n", float(elapsed_seconds)/1000.);
		}
	
		
	
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
		map<string, vector<merged_interval*> > G 	= load_intervals(interval_file, pad); //load the forward and reverse strand intervals, merge accordingly...
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
		MPI::Finalize();
		return 1;
	}
	delete P;
	MPI::Finalize();
	
	return 0;
}
