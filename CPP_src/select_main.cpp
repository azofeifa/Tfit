#include "load.h"
#include "select_main.h"
#include "MPI_comm.h"
#include "template_matching.h"
#include <omp.h>
#include "model_main.h"
#include "model_selection.h"
int select_run(params * P, int rank, int nprocs, int job_ID, Log_File * LG){
	int verbose 	= stoi(P->p["-v"]);

	LG->write("\ninitializing select module..............................done\n\n",verbose);

	int threads 	= omp_get_max_threads();//number of OpenMP threads that are available for use
	string job_name = P->p["-N"];
	
	
	string forward_bedgraph 		= P->p["-i"]; //forward strand bedgraph file
	string reverse_bedgraph 		= P->p["-j"]; //reverse strand bedgraph file
	string out_file_dir 			= P->p["-o"];//out file directory
	string query_file 				= P->p["-q"];
	string tss_file 				= P->p["-tss"];
	
	map<string, int> chrom_to_ID;
	map<int, string> ID_to_chrom;

	LG->write("loading bedgraph files..................................",verbose);
	vector<segment *> 	segments 	= load::load_bedgraphs_total(forward_bedgraph, reverse_bedgraph, stoi(P->p["-br"]), stof(P->p["-ns"]), P->p["-chr"], chrom_to_ID, ID_to_chrom );
	vector<segment *> all_segments  	= segments;
	LG->write("done\n",verbose);

	LG->write("slicing segments........................................",verbose);
	segments 						= MPI_comm::slice_segments(segments, rank, nprocs);	
	LG->write("done\n",verbose);
	
	//====================================================================
	//(2a) get noise intervals, bct:LLR < 0.2, get how many?
	//this will pretty much follow the bidir module pipeline
	LG->write("global template matching for non-divergent signal.......",verbose);
	noise_global_template_matching(segments, stod(P->p["-ns"]) );
	LG->write("done\n",verbose);

	//(2b) gather and write out
	LG->write("gather non-divergent signal.............................",verbose);
	int total =  MPI_comm::gather_all_bidir_predicitions(all_segments, 
			segments , rank, nprocs, out_file_dir, job_name, job_ID,P,1);
	LG->write("done\n",verbose);

	LG->write("clearing allocated segment memory.......................",verbose);	
	load::clear_segments(all_segments);
	LG->write("done\n",verbose);
	//====================================================================
	//(3) run model across noise intervals
	//this will pretty much follow the model module pipeline
	//out_dir+ job_name+ "-" + to_string(job_ID)+ "_prelim_bidir_hits.bed"
	P->p["-k"] 	= out_file_dir+ job_name+ "-" + to_string(job_ID)+ "_non_div_txn.bed";
	P->p["-minK"] 	= "1",P->p["-maxK"] 	= "1";
	model_run(P, rank, nprocs, 1,  job_ID, LG);


	//====================================================================
	//(4) load the noise_K_models and the query K_models
	//query 
	string noise_K_models_file 	= out_file_dir +  P->p["-N"] + "-" + to_string(job_ID)+  "_non_div_txn_K_models_MLE.tsv";
	string noise_bidir_ms_file 	= P->p["-o"]+  P->p["-N"] + "-" + to_string(job_ID)+  "_non_div_txn_divergent_classifications.bed";
	double AUC = 0, TP = 0, FP = 0, TP_at_fp=0,FP_at_tp=0,optimal_penality=1;
	if (rank==0){
		LG->write("loading K_model out files...............................",verbose);
		vector<segment_fits *> query_fits 		= load::load_K_models_out(query_file);
		vector<segment_fits *> noise_fits 		= load::load_K_models_out(noise_K_models_file);
		LG->write("done\n",verbose);
		LG->write("labeling TSS overlap (true positive assumption).........",verbose);
		vector<segment_fits *> query_fits_tss 	= load::label_tss(tss_file, query_fits);
		LG->write("done\n",verbose);
		LG->write("calculating ROC curve and optimizing penalized BIC......",verbose);
		ROC(noise_fits, query_fits, AUC, TP, FP, TP_at_fp, FP_at_tp, optimal_penality);
		LG->write("done\n",verbose);

	}

	//====================================================================
	//(5) need to label K_models as to whether they are 
	// over transcription start sites (these are true positives)
	// compute ROC at get BIC penalty for the TP or FP rate

	
	if (rank == 0){
		LG->write("deleting noise model files..............................",verbose);
		remove( P->p["-k"].c_str() ) ;
		remove( noise_K_models_file.c_str() ) ;
		remove( noise_bidir_ms_file.c_str() );
		LG->write("done\n\n",verbose);
		LG->write("AUC              : " + to_string(AUC) + "\n", verbose);
		LG->write("Optimal TP       : " + to_string(TP) + "\n", verbose);
		LG->write("Optimal FP       : " + to_string(FP) + "\n", verbose);
		LG->write("Optimal Penalty  : " + to_string(optimal_penality) + "\n", verbose);
	}

	LG->write("\nexiting select module...................................done\n\n",verbose);
	//====================================================================
	// (6) finally output bidirectional output file 

	P->p["-ms_pen"] 	= to_string(optimal_penality);

	MPI_comm::wait_on_root(rank, nprocs);

	
}
