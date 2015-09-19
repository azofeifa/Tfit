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
#include "omp.h"
using namespace std;

class timer{
public:
	clock_t t;		
	chrono::time_point<chrono::system_clock> start, end;
	int WT;
	string HEADER;
	void start_time(int rank, string header){
		if (rank==0){
			string white_space 	= "";
			if (header.size() > WT){
				WT 				= header.size();
			}
			for (int i = 0; i < (WT-header.size() ); i++ ){
				white_space+= " ";
			}
			HEADER=header+white_space;
			start = chrono::system_clock::now();
			t = clock();
		}
	}
	string get_time(int rank){
		string out 	= "";
			
		if (rank==0){
			end = chrono::system_clock::now();
			int elapsed_seconds = std::chrono::duration_cast<std::chrono::milliseconds >
			                             (end-start).count();
			t = clock() - t;
			if (rank ==0 ){
				out+=HEADER;
				out+="CPU: " + to_string(t / double(CLOCKS_PER_SEC));
				out+="WT:  " + to_string(elapsed_seconds/1000.);
				printf("%s",HEADER.c_str() );
				printf("CPU: %.5g,",t / double(CLOCKS_PER_SEC)  );
				printf("WT: %.5g\n", (elapsed_seconds/1000.));
			}
		}
		return out;

	}
	timer(int wt){
		WT = wt;
	}
};


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
		int nprocs		= MPI::COMM_WORLD.Get_size();
		int rank 		= MPI::COMM_WORLD.Get_rank();
	    int threads  	= omp_get_max_threads();
		int verbose 	= stoi(P->p4["-v"]);
		string log_out 	= P->p4["-log_out"] + "tmp_log_file_" + to_string(rank) + ".log"  ;
		ofstream 	FHW;
		FHW.open(log_out);

		if (verbose and rank==0){//show current user parameters...
			P->display(nprocs,threads);
		}
		FHW<<"#Temp Log File for mpi process: " + to_string(rank) + "\n";
		FHW<<P->get_header(4);


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
		double window 				= stod(P->p4["-window_res"]);
		double ct 					= stod(P->p4["-bct"]);
		double density 				= stod(P->p4["-density"]) / scale;
		int opt_res 				= stod(P->p4["-opt_res"]);
		int np 						= stoi(P->p4["-np"]);
		char processor_name[MPI_MAX_PROCESSOR_NAME];
		int namelen;  
		MPI_Get_processor_name(processor_name, &namelen); 
		string node_name 			= processor_name;
		string spec_chrom 			= P->p4["-chr"];
		timer T(80);
		timer TF(80);
		
		TF.start_time(rank, "Final Time:");
		
		T.start_time(rank, "loading BG files:");
		FHW<<"(main) loaded begraph files...";
		vector<segment*> segments 	= load_bedgraphs_total(forward_bedgraph, 
			reverse_bedgraph, BINS, scale, spec_chrom);
		FHW<<"done\n";
		T.get_time(rank);
		if (segments.empty()){
			printf("segments not populated, exiting...\n");
			MPI::Finalize();
			return 1;
		}
		vector<segment*> all_segments  	= segments;
		
		segments 						= slice_segments(segments, rank, nprocs);
		FHW<<"(main) sliced segments: " + to_string(int(segments.size()))+ " on this process\n";
		FHW.flush();
		if (not optimize_directory.empty() and rank==0 and not segments.empty() ){
			map<string, interval_tree *> I = load_bidir_bed_files(optimize_directory, 
				spec_chrom);

			if (I.empty()){
				printf("exiting...\n");
				MPI::Finalize();
				return 1;
			}
			free_segments(segments);
			MPI::Finalize();
			return 0;			
		}else if(optimize_directory.empty() and not segments.empty() ){
			T.start_time(rank, "running template matching:");	
			run_global_template_matching(segments, out_file_dir, window, 
				density,scale,ct, np,0. ,0, FHW );	
			T.get_time(rank);
			FHW<<"ran template matching algorithm\n";
			
		}
		map<string , vector<vector<double> > > G;
		T.start_time(rank, "(MPI) gathering bidir predictions:");	
		FHW<<"(main) gathering all bidir predictions...";
		if (P->p4["-show_seeds"] == "1"){
			G = gather_all_bidir_predicitions(all_segments, 
				segments , rank, nprocs, out_file_dir);
		}else{
			G = gather_all_bidir_predicitions(all_segments, 
				segments , rank, nprocs, "");
		}
		FHW<<"done\n";
		FHW.flush();
		
		T.get_time(rank);

		if (P->p4["-MLE"] == "1"){
			vector<segment *> bidir_segments;
			if (not G.empty()  ){
				T.start_time(rank, "loading BG files:");
				FHW<<"(main) loading bidir prediction files...";
				FHW.flush();

				bidir_segments 	= bidir_to_segment( G, 
					forward_bedgraph,reverse_bedgraph, stoi(P->p4["-pad"]),P->p4["-chr"]   );
				T.get_time(rank);
				FHW<<"done\n";

				FHW.flush();

			}
			vector<simple_c> fits;
			if (not bidir_segments.empty()){
				FHW<<"(main) bining bidir segments...";
				BIN(bidir_segments, stod(P->p4["-br"]), stod(P->p4["-ns"]),true );
				FHW<<"done\n";
				FHW.flush();
				T.start_time(0, "MLE fit on " + node_name + ", going to process " + to_string(int(bidir_segments.size())) + " segments: ");
				fits 			= run_model_accross_segments_to_simple_c(bidir_segments, P,FHW);
				T.get_time(0);
			}
			FHW<<"(main) (MPI) gathering MLE results...";
				
			T.start_time(rank, "(MPI) gathering MLE results:");
			map<string, map<int, vector<rsimple_c> > > rcG 	= gather_all_simple_c_fits(bidir_segments, fits, rank, nprocs);
			FHW<<"done";
			FHW.flush();
			T.get_time(rank);
			
			vector<segment *> FSI;
			map<int, string> IDS;	
					
			if (rank==0 and not rcG.empty() ){//perform and optimize model selection based on number of bidir counts
				T.start_time(rank, "opt model selection:");
				vector<final_model_output> 	A  				= optimize_model_selection_bidirs(rcG, P, FHW);
				T.get_time(rank);
				T.start_time(rank, "writing out bidir model selection:");
				write_out_MLE_model_info(A, P);
				T.get_time(rank);
				

				if (P->p4["-elon"] == "1" and rank==0){
					//want load the intervals of "interest"
					T.start_time(rank, "combinding FSI and bidir intervals:");
					FSI 		= load_intervals_of_interest(P->p4["-f"], IDS ,0 );
					//now we want to insert final_model_output data into FSI...	
					FSI 		= combind_bidir_fits_with_intervals_of_interest( A,  FSI );		
					T.get_time(rank);			
				}
			}
			fits.clear();
			if (P->p4["-elon"] == "1"){
				FHW<<"(main) (MPI) sending out elongation assignments...";
				T.start_time(rank, "(MPI) sending out elongation assignments:");
				map<string, vector<segment *> > GG 	= send_out_elongation_assignments(FSI, rank, nprocs);
				T.get_time(rank);
				FHW<<"done\n";
				FHW.flush();
				vector<segment*> integrated_segments;
				if (not GG.empty()){
					T.start_time(rank, "loading BG and integrating segments:");
					integrated_segments= insert_bedgraph_to_segment(GG, forward_bedgraph ,reverse_bedgraph,rank);
					printf("Rank: %d, going to move elongation support on %d \n",rank, int(integrated_segments.size() )  );
					BIN(integrated_segments, stod(P->p4["-br"]), stod(P->p4["-ns"]),true);
					T.get_time(rank);
					T.start_time(0, "shift elon. rank " + to_string(rank) + ",on " + to_string(int(integrated_segments.size())) + " segments:");
					fits 	= move_elongation_support(integrated_segments, P);
					T.get_time(0);
				}
				T.start_time(rank, "(MPI) gathering elongation results:");
				map<string, map<int, vector<rsimple_c> > > rcG 	= gather_all_simple_c_fits(integrated_segments, fits, rank, nprocs);
				vector<final_model_output> A 					= convert_to_final_model_output(rcG, P);
				//convert to final_model_output
				if (rank==0){
					T.get_time(rank);
					T.start_time(rank, "Writing Out Model Fits:");
					write_gtf_file_model_fits(A, P);
					write_config_file_model_fits(A, IDS, P);
					T.get_time(rank);
				}
				
			}
		}
		if (rank==0){
			collect_all_tmp_files(P->p4["-log_out"], nprocs);
		}
		

		if (not segments.empty()){
			free_segments(segments);
		}
		MPI::Finalize();
		TF.get_time(rank);

		return 0;
    
	}else if (P->module=="SINGLE"){
		int nprocs = MPI::COMM_WORLD.Get_size();
		int rank = MPI::COMM_WORLD.Get_rank();
		int verbose 	= stoi(P->p3["-v"]);
	   	int threads  	= omp_get_num_threads();

		if (verbose and rank==0){//show current user parameters...
			P->display(nprocs, threads);
		}	
		string bed_graph_file 	= P->p5["-i"];
		string interval_file 	= P->p5["-j"];
		string out_file_dir 		= P->p5["-o"];
		int BINS 					= stoi(P->p5["-br"]);
		double scale 				= stod(P->p5["-ns"]);
		double ct 					= stod(P->p5["-bct"]);
		double opt_res 				= stod(P->p5["-opt_res"]);
		int np 						= stoi(P->p5["-np"]);
		
		string spec_chrom 			= P->p5["-chr"];
		bool run_template 			= bool(stoi(P->p5["-template"]));
		vector<segment *> FSI;
		map<string , vector<vector<double> > > G;
		map<string, vector<segment *> > GG;
		timer T(50);
				
		map<int, string> IDS;
		if (rank==0){
			T.start_time(rank, "Loading/Converting intervals of interest:");
			FSI 							= load_intervals_of_interest(interval_file, IDS, stoi(P->p5["-pad"]) );
			if (not G.empty()){
				vector<final_model_output> 	A 	= convert_bidir_segs_to_final_model(G);
				combind_bidir_fits_with_intervals_of_interest( A,  FSI );		
			}
			T.get_time(rank);					
		}
		T.start_time(rank, "(MPI) sending out elongation assignments:");
		GG 	= send_out_single_fit_assignments(FSI, rank, nprocs);
		T.get_time(rank);
		vector<segment*> integrated_segments;
		vector<single_simple_c> single_fits;
		if (not GG.empty()){
			T.start_time(rank, "loading BG and integrating segments:");
			integrated_segments= insert_bedgraph_to_segment_single(GG, bed_graph_file,rank);
			BIN(integrated_segments, stod(P->p5["-br"]), stod(P->p5["-ns"]),true);
			T.get_time(rank);

			//now running model
			T.start_time(rank, "running model accross segments:");
			single_fits = run_single_model_across_segments(integrated_segments, P);
			T.get_time(rank);
		}
		T.start_time(rank, "(MPI) Gathering Single Fit Info:");
		vector<single_simple_c> all_fits 	= gather_all_simple_c(single_fits, rank, nprocs);
		T.get_time(rank);
		
		if (rank == 0 and not all_fits.empty()){
			write_out_single_simple_c(all_fits, IDS , P );
		}



					




		
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
