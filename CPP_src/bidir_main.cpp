#include "bidir_main.h"
#include "template_matching.h"
#include "density_profiler.h"
#include "MPI_comm.h"
#include <omp.h>
#include "model_main.h"
#include "select_main.h"
#include "error_stdo_logging.h"
using namespace std;
int bidir_run(params * P, int rank, int nprocs, int job_ID, Log_File * LG){
	int verbose 	= stoi(P->p["-v"]);
	P->p["-merge"] 	= "1";
	LG->write("\ninitializing bidir module...............................done\n", verbose);
	int threads 	= omp_get_max_threads();//number of OpenMP threads that are available for use
	
	//===========================================================================
	//get job_ID and open file handle for log files
	string job_name = P->p["-N"];
	
	//===========================================================================
	//input files and output directories
	string forward_bedgraph 	= P->p["-i"]; //forward strand bedgraph file
	string reverse_bedgraph 	= P->p["-j"]; //reverse strand bedgraph file
	string out_file_dir 		= P->p["-o"] ;//out file directory
	//(2a) read in bedgraph files 
	map<string, int> chrom_to_ID;
	map<int, string> ID_to_chrom;
	LG->write("loading bedgraph files..................................", verbose);
	vector<segment *> 	segments 	= load::load_bedgraphs_total(forward_bedgraph, reverse_bedgraph, 
		stoi(P->p["-br"]), stof(P->p["-ns"]), P->p["-chr"], chrom_to_ID, ID_to_chrom );
	LG->write("done\n", verbose);
	//(2b) so segments is indexed by inidividual chromosomes, want to broadcast 
	//to sub-processes and have each MPI call run on a subset of segments
	vector<segment*> all_segments  	= segments;
	LG->write("slicing segments........................................", verbose);
	segments 						= MPI_comm::slice_segments(segments, rank, nprocs);	
	LG->write("done\n", verbose);
	//===========================================================================
	//(3a) now going to run the template matching algorithm based on pseudo-
	//moment estimator and compute BIC ratio (basically penalized LLR)

	LG->write("running moment estimator algorithm......................", verbose);
	run_global_template_matching(segments, out_file_dir, 4, 
			0.,stod(P->p["-ns"]),stod(P->p["-bct"]), threads,0. ,0 );	
	//(3b) now need to send out, gather and write bidirectional intervals 
	LG->write("done\n", verbose);
	LG->write("scattering predictions to other MPI processes...........", verbose);
	int total =  MPI_comm::gather_all_bidir_predicitions(all_segments, 
			segments , rank, nprocs, out_file_dir, job_name, job_ID,P,0);
	LG->write("done\n", verbose);
	if (rank==0){
		LG->write("\nThere were " +to_string(total) + " prelimary bidirectional predictions\n\n", verbose);
	}
	
	//===========================================================================
	//this should conclude it all
	LG->write("clearing allocated segment memory.......................", verbose);	
	load::clear_segments(all_segments);
	LG->write("done\n", verbose);
	//===========================================================================
	//(4) if MLE option was provided than need to run the model_main::run()
	//
	if (stoi(P->p["-MLE"])){
		P->p["-k"] 	= P->p["-o"]+ job_name+ "-" + to_string(job_ID)+ "_prelim_bidir_hits.bed";
		model_run(P, rank, nprocs,0, job_ID, LG);
		if (stoi(P->p["-select"])){
			string bidir_ms_file 	= P->p["-o"]+  P->p["-N"] + "-" + to_string(job_ID)+  "_divergent_classifications.bed";
			string file_name 		= P->p["-o"]+  P->p["-N"] + "-" + to_string(job_ID)+  "_K_models_MLE.tsv";
			remove( bidir_ms_file.c_str() );
			//query out_dir +  P->p["-N"] + "-" + to_string(job_ID)+  "_K_models_MLE.tsv"
			P->p["-q"] 	= P->p["-o"] +  P->p["-N"] + "-" + to_string(job_ID)+  "_K_models_MLE.tsv";
			select_run(P, rank, nprocs, job_ID, LG);
			LG->write("loading results (MLE)...................................",verbose);
			vector<segment_fits *> fits 		= load::load_K_models_out(file_name);
			LG->write("done\n",verbose);		
			
			LG->write("writing out results (model selection)...................",verbose);
			load::write_out_bidirectionals_ms_pen(fits, P,job_ID, 0 );
			LG->write("done\n",verbose);
	

		}
	}
	LG->write("exiting bidir module....................................done\n\n", verbose);
	return 1;
}
