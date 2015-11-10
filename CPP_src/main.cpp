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
#include "bootstrap.h"
#include "density_profiler.h"
#include <omp.h>
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
		string job_name = P->p4["-N"];
		int job_ID 		=  get_job_ID(P->p4["-log_out"], job_name, rank, nprocs);
		string log_out 	= P->p4["-log_out"] + "tmp_" + job_name+ "-"+ to_string(job_ID)+ "_" + to_string(rank) + ".log"  ;
		ofstream 	FHW;
		FHW.open(log_out);

		if (verbose and rank==0){//show current user parameters...
			P->display(nprocs,threads);
		}
		FHW<<"#Temp Log File for mpi process: " + to_string(rank) + "\n";
		FHW<<P->get_header(4);


		//load bed graph files into map<string, double **>;
		string forward_bedgraph 	= P->p4["-i"];
		string reverse_bedgraph 	= P->p4["-j"];
		string noise_bed_file 		= P->p4["-nf"];
		string out_file_dir 		= P->p4["-o"] ;
		
		int BINS 					= stoi(P->p4["-br"]);
		double scale 				= stod(P->p4["-ns"]);
		double window 				= stod(P->p4["-window_res"]);
		double ct 					= stod(P->p4["-bct"]);
		int opt_res 				= stod(P->p4["-opt_res"]);
		int np 						= stoi(P->p4["-np"]);
		char processor_name[MPI_MAX_PROCESSOR_NAME];
		int namelen;  
		MPI_Get_processor_name(processor_name, &namelen); 
		string node_name 			= processor_name;
		string spec_chrom 			= P->p4["-chr"];
		timer T(80);
		timer TF(80);
		
		double mean = 0;
		double var 	= 0;

		if (rank==0){
		  	//get_noise_mean_var(noise_bed_file, forward_bedgraph, &mean, &var);
		  	mean 	=  get_table_mean_var(noise_bed_file, forward_bedgraph, window, stof(P->p4["-br"]), scale );
		}
		double density 		= send_density_val((mean  )*scale , rank, nprocs );
		string NOISE_OUT 	= "(main) estimated noise density: " + to_string(density/scale);
		if (rank==0){
			printf("%s\n", NOISE_OUT.c_str() );
			FHW<<NOISE_OUT<<endl;
		}
		TF.start_time(rank, "Final Time:");
		
		T.start_time(rank, "loading BG files:");
		FHW<<"(main) loaded begraph files...";
		map<string, int> chrom_to_ID;
		map<int, string> ID_to_chrom;
		vector<segment*> segments 	= load_bedgraphs_total(forward_bedgraph, 
			reverse_bedgraph, BINS, scale, spec_chrom, chrom_to_ID, ID_to_chrom );
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
		if(  not segments.empty() ){
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
				segments , rank, nprocs, out_file_dir, job_name, job_ID,P,FHW);
		}else{
			G = gather_all_bidir_predicitions(all_segments, 
				segments , rank, nprocs, "", job_name, job_ID,P,FHW);
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
					forward_bedgraph,reverse_bedgraph, stoi(P->p4["-pad"]),P->p4["-chr"],chrom_to_ID   );
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
				FHW<< "MLE fit on " + node_name + ", going to process " + to_string(int(bidir_segments.size())) + " segments"<<endl;
				FHW.flush();

				T.start_time(0, "MLE fit on " + node_name + ", going to process " + to_string(int(bidir_segments.size())) + " segments: ");
				fits 			= run_model_accross_segments_to_simple_c(bidir_segments, P,FHW);
				T.get_time(0);
			}
			FHW<<"(main) (MPI) gathering MLE results...";
				
			T.start_time(rank, "(MPI) gathering MLE results:");
			map<string, map<int, vector<rsimple_c> > > rcG 	= gather_all_simple_c_fits(bidir_segments, 
				fits, rank, nprocs,ID_to_chrom);
			FHW<<"done";
			FHW.flush();
			T.get_time(rank);
			
			vector<segment *> FSI;
			map<int, string> IDS;	
					
			if (rank==0 and not rcG.empty() ){//perform and optimize model selection based on number of bidir counts
				T.start_time(rank, "opt model selection:");
				vector<final_model_output> 	A  				= optimize_model_selection_bidirs(rcG, P, FHW,ID_to_chrom);
				T.get_time(rank);
				T.start_time(rank, "writing out bidir model selection:");
				write_out_MLE_model_info(A, P, job_name, job_ID);
				T.get_time(rank);
				

			}
			fits.clear();
			
		}
		if (rank==0){
			collect_all_tmp_files(P->p4["-log_out"], job_name, nprocs, job_ID);
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
	   	string job_name = P->p5["-N"];
	   	cout<<P->p5["-log_out"]<<endl;
		int job_ID 		=  get_job_ID(P->p5["-log_out"], job_name, rank, nprocs);
		string log_out 	= P->p5["-log_out"] + "tmp_" + job_name+ "-"+ to_string(job_ID)+ "_" + to_string(rank) + ".log"  ;
		ofstream 	FHW;
		FHW.open(log_out);
		FHW<<"#Temp Log File for mpi process: " + to_string(rank) + "\n";
		FHW<<P->get_header(5);
		FHW.flush();
		if (verbose and rank==0){//show current user parameters...
			P->display(nprocs, threads);
		}	
		string bed_graph_file 	= P->p5["-i"];
		string interval_file 	= P->p5["-j"];
		string out_file_dir 		= P->p5["-o"];
		int BINS 					= stoi(P->p5["-br"]);
		double scale 				= stod(P->p5["-ns"]);
		
		string spec_chrom 			= P->p5["-chr"];
		vector<segment *> FSI;
		map<string , vector<vector<double> > > G;
		map<string, vector<segment *> > GG;
		timer T(50);
				
		map<int, string> IDS;
		if (rank==0){
			T.start_time(rank, "Loading/Converting intervals of interest:");
			FSI 							= load_intervals_of_interest(interval_file, 
				IDS, stoi(P->p5["-pad"]), spec_chrom );
			FHW<<"Loading/Converting intervals of interest"<<endl;
			FHW.flush();

			if (stoi(P->p5["-merge"]) ){
				FSI 						= merge_intervals_of_interest(FSI);
			}
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
			FHW<<"binned and loaded intervals of integrated segments"<<endl;
			FHW.flush();

			//now running model
			T.start_time(rank, "running model accross segments:");
			single_fits = run_single_model_across_segments(integrated_segments, P, FHW);
			T.get_time(rank);
		}
		T.start_time(rank, "(MPI) Gathering Single Fit Info:");
		vector<single_simple_c> all_fits 	= gather_all_simple_c(single_fits, rank, nprocs);
		T.get_time(rank);
		if (rank==0){
			write_out_single_simple_c_marks(all_fits, IDS , P );
			collect_all_tmp_files(P->p5["-log_out"], job_name, nprocs, job_ID);
		}
		

		
		
	}else if (P->module=="MODEL"){
		int nprocs = MPI::COMM_WORLD.Get_size();
		int rank = MPI::COMM_WORLD.Get_rank();
		int verbose 	= stoi(P->p["-v"]);
	   	int threads  	= omp_get_max_threads();
		char processor_name[MPI_MAX_PROCESSOR_NAME];
		int namelen;  
		MPI_Get_processor_name(processor_name, &namelen); 
		string job_name = P->p["-N"];

		int job_ID 		=  get_job_ID(P->p["-log_out"], job_name, rank, nprocs);
		string log_out 	= P->p["-log_out"] + "tmp_" + job_name+ "-"+ to_string(job_ID)+ "_" + to_string(rank) + ".log"  ;
		ofstream 	FHW;
		FHW.open(log_out);
		timer TF(50);
		
		TF.start_time(rank, "Final Time:");
		
		
		if (verbose and rank==0){//show current user parameters...
			P->display(nprocs, threads);
		}	
		if (rank==0){
			FHW<<P->get_header(0);
		}
		FHW.flush();
		
		string forward_bed_graph_file 	= P->p["-i"];
		string reverse_bed_graph_file 	= P->p["-j"];

		string interval_file 			= P->p["-k"];
		string out_file_dir 			= P->p["-o"];
		string spec_chrom 			= P->p["-chr"];
		bool run_template 			= bool(stoi(P->p["-template"]));
		vector<segment *> FSI;
		map<string , vector<vector<double> > > G;
		map<string, vector<segment *> > GG;
		timer T(50);
		

		map<int, string> IDS;
		if (rank==0){
			FHW<<"(main) Loading/Converting intervals of interest"<<endl;
			T.start_time(rank, "Loading/Converting intervals of interest:");
			FSI 							= load_intervals_of_interest(interval_file, IDS, stoi(P->p["-pad"]), spec_chrom );
			T.get_time(rank);					
			FHW.flush();
		}
		T.start_time(rank, "(MPI) sending out segment assignments:");
		FHW<<"(main) waiting to receive elongation assignments"<<endl;
		FHW.flush();
		GG 	= send_out_single_fit_assignments(FSI, rank, nprocs);
		vector<segment*> integrated_segments;
		FHW<<"(main) received elongation assignments"<<endl;
		FHW.flush();

		T.get_time(rank);
		FHW<<"(main) Loading bedgraph files into intervals of interest: ";
		integrated_segments= insert_bedgraph_to_segment_joint(GG, 
			forward_bed_graph_file, reverse_bed_graph_file, rank);
		FHW<<to_string(int(integrated_segments.size()))<<endl;
		FHW.flush();
		FHW<<"(main) binning integrated segments: ";
		FHW.flush();
		BIN(integrated_segments, stod(P->p["-br"]), stod(P->p["-ns"]),true);
		FHW<<"DONE"<<endl;
		FHW.flush();
			
		double window 				= stod(P->p["-window_res"]);
		double ct 					= stod(P->p["-bct"]);
		double scale 				= stod(P->p["-ns"]);

		T.start_time(rank, "Running Template Matching on individual segments:");
		FHW.flush();
		run_global_template_matching(integrated_segments, out_file_dir, window, 
				0.1,scale,ct, 64,0. ,0, FHW );	
		T.get_time(rank);
		T.start_time(rank, "Running Model on individual segments:");
		FHW<<"(main) About to run model"<<endl;
		FHW.flush();
		vector<map<int, vector<simple_c_free_mode> >> FITS 		= run_model_across_free_mode(integrated_segments,
		 P,FHW);
		T.get_time(rank);
		FHW<<"(main) Gathering all simple c free model"<<endl;		
		FHW.flush();
		map<int, map<int, vector<simple_c_free_mode>  > > GGG 	= gather_all_simple_c_free_mode(FITS, rank, nprocs);
		if (rank==0){//write_out_to_MLE
			write_out_models_from_free_mode(GGG, P, job_ID, IDS);
		}
		FHW<<"(main) gathered all simple c free model"<<endl;		
		FHW.flush();
		
		if (rank==0){
			collect_all_tmp_files(P->p["-log_out"], job_name, nprocs, job_ID);
		}
		TF.get_time(rank);		
	}
	else if (P->module=="BOOTSTRAP"){
		int nprocs		= MPI::COMM_WORLD.Get_size();
		int rank 		= MPI::COMM_WORLD.Get_rank();
	    int threads  	= omp_get_max_threads();
		int verbose 	= stoi(P->p6["-v"]);
		string job_name = P->p6["-N"];
		fill_in_bidir_boostrap(P);
		int job_ID 		= get_job_ID(P->p6["-log_out"], job_name, rank, nprocs);

		map<string, int> chrom_to_ID;
		map<int, string> ID_to_chrom;
		
		string log_out 	= P->p6["-log_out"] + "tmp_" + job_name+ "-"+ to_string(job_ID)+ "_" + to_string(rank) + ".log"  ;
		ofstream 	FHW;
		FHW.open(log_out);
		FHW<<"#Temp Log File for mpi process: " + to_string(rank) + "\n";
		FHW<<P->get_header(6);
		if (verbose and rank==0){//show current user parameters...
			P->display(nprocs,threads);
		}
		vector<vector<int> > start_stops;
		if (rank==0){
			printf("getting line start and stops\n");
			start_stops 	=  get_line_start_stops(P, nprocs);
		}
		if (rank==0){
			printf("(MPI) sending out line start and stops\n");
		}
		vector<int> st_sp 						= send_out_merged_start_stops(start_stops,   rank,   nprocs);
		if (rank==0){
			printf("loading bidirectional predictions\n");
		}
		map<string, vector<segment *> > GG		= load_bidir_predictions( P, 
			st_sp,chrom_to_ID,ID_to_chrom) ;


		if (rank==0){
			printf("inserting bedgraph data\n");
		}
		vector<segment *> integrated_segments= insert_bedgraph_to_segment_joint(GG, 
			P->p6["-j"], P->p6["-k"], rank);
		if (rank==0){
			printf("binning\n");
		}
		BIN(integrated_segments, stod(P->p6["-br"]), stod(P->p6["-ns"]),true);
		if (rank==0){
			printf("running bootstrap\n");
		}
		run_bootstrap_across(integrated_segments, P, FHW);
		vector<boostrap_struct> bs 	= collect_bootstrap(integrated_segments,  rank,  nprocs, chrom_to_ID);
		if (rank == 0){
			write_bootstrap( bs,  ID_to_chrom,  P, job_ID);
		}
		if (rank==0){
			collect_all_tmp_files(P->p6["-log_out"], job_name, nprocs, job_ID);
		}
	
		FHW.close();
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
