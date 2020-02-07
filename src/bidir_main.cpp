#include "bidir_main.h"
#include "template_matching.h"
#include "density_profiler.h"
#include "MPI_comm.h"
#include <mpi.h>
#include <omp.h>
#include "model_main.h"
#include "select_main.h"
#include "error_stdo_logging.h"
#include "FDR.h"
#include "BIC.h"
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
	string joint_bedgraph 		= P->p["-ij"]; //joint forward and reverse strand bedgraph file
	string tss_file 		= P->p["-tss"];
	string out_file_dir 		= P->p["-o"] ;//out file directory
	//===========================================================================
	//template searching parameters
	double sigma, lambda, foot_print, pi, w;
	double ns 	= stod(P->p["-ns"]);
	sigma 		= stod(P->p["-sigma"]), lambda= stod(P->p["-lambda"]);
	foot_print 	= stod(P->p["-foot_print"]), pi= stod(P->p["-pi"]), w= stod(P->p["-w"]);
	//(2a) read in bedgraph files 
	map<string, int> chrom_to_ID;
	map<int, string> ID_to_chrom;
       
	vector<double> parameters 	= {sigma, lambda, foot_print,pi, w};
	if (not tss_file.empty() and rank == 0){
		vector<segment *> FSI;
		LG->write("loading TSS intervals...................................",verbose);
		map<int, string> IDS;
		vector<segment *> tss_intervals 	= load::load_intervals_of_interest(tss_file, 
			IDS, P, 1 );
		map<string, vector<segment *>> GG 	= MPI_comm::convert_segment_vector(tss_intervals);
		LG->write("done\n", verbose);
		
		LG->write("inserting coverage data.................................",verbose);
		vector<segment*> integrated_segments= load::insert_bedgraph_to_segment_joint(GG, 
			forward_bedgraph, reverse_bedgraph, joint_bedgraph, rank);
		LG->write("done\n", verbose);

		LG->write("Binning/Normalizing TSS intervals.......................",verbose);
		load::BIN(integrated_segments, stod(P->p["-br"]), stod(P->p["-ns"]),true);	
		LG->write("done\n", verbose);
		
		LG->write("Computing Average Model.................................",verbose);
		parameters  	= compute_average_model(integrated_segments, P);
		//these values are on the genome scale (non normalize to ns)
		LG->write("done\n", verbose);
		LG->write("\nAverage Model Parameters\n", verbose);
		LG->write("-sigma      : " + to_string(parameters[0])+ "\n", verbose);
		LG->write("-lambda     : " + to_string(parameters[1])+ "\n", verbose);
		LG->write("-foot_print : " + to_string(parameters[2])+ "\n", verbose);
		LG->write("-pi         : " + to_string(parameters[3])+ "\n", verbose);
		LG->write("-w          : " + to_string(parameters[4])+ "\n\n", verbose);
	}
	parameters 				= MPI_comm::send_out_parameters( parameters, rank, nprocs);
	P->p["-sigma"] 	       = to_string(parameters[0]);
	P->p["-lambda"]        = to_string(parameters[1]);
	P->p["-foot_print"]    = to_string(parameters[2]);
	P->p["-pi"] 	       = to_string(parameters[3]);
	P->p["-w"] 	       = to_string(parameters[4]);

	LG->write("loading bedgraph files..................................", verbose);
	vector<segment *> 	segments 	= load::load_bedgraphs_total(forward_bedgraph, reverse_bedgraph, joint_bedgraph,
		stoi(P->p["-br"]), stof(P->p["-ns"]), P->p["-chr"], chrom_to_ID, ID_to_chrom );

	if (segments.empty()){
		printf("exiting...\n");
		return 1;
	}
	LG->write("done\n", verbose);




	slice_ratio SC;
	if (stoi(P->p["-FDR"] ) ){
	  LG->write("getting likelihood score distribution...................", verbose);
	  SC                      = get_slice(segments, pow(10,6) , pow(10,4) ,P  );
	  LG->write("done\n\n", verbose);
	  if (not SC.converged){
	    LG->write("converged            : False (restoring default values)\n"  ,verbose );
	  }else{
	    LG->write("converged            : True\n" ,verbose );
	  }
	  LG->write("score mean           : "+to_string(SC.mean) + "\n" ,verbose );
	  LG->write("standard Deviation   : "+to_string(SC.std ) + "\n" ,verbose );
	  LG->write("h                    : "+to_string(SC.w ) + "\n" ,verbose );
	  LG->write("threshold            : "+to_string(SC.threshold) + "\n\n" ,verbose );
	}
	else{
	  SC.mean = 0.6, SC.std = 0.001; //this dependent on -w 0.9 !!!
	  SC.set_2(stod(P->p["-bct"]));
	}

	//================================================
	//computing this SC via multiple times of the EM
	//run, should probably write something to combind
	//them all; i.e. MPI

	//=================================================
	//(2b) so segments is indexed by inidividual chromosomes, want to broadcast 
	//to sub-processes and have each MPI call run on a subset of segments
	vector<segment*> all_segments  	= segments;
	LG->write("slicing segments........................................", verbose);
	segments 						= MPI_comm::slice_segments(segments, rank, nprocs);	
	LG->write("done\n", verbose);

	//===========================================================================
	//(3a) now going to run the template matching algorithm based on pseudo-
	//moment estimator and compute BIC ratio (basically penalized LLR)
	LG->write("running template matching algorithm.....................", verbose);
	double threshold 	= run_global_template_matching(segments, out_file_dir, P, SC);	
	//(3b) now need to send out, gather and write bidirectional intervals 
	LG->write("done\n", verbose);
	


	LG->write("scattering predictions to other MPI processes...........", verbose);
	int total =  MPI_comm::gather_all_bidir_predicitions(all_segments, 
							     segments , rank, nprocs, out_file_dir, job_name, job_ID,P,0);
	MPI_Barrier(MPI_COMM_WORLD); //make sure everybody is caught up!

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
		P->p["-k"] 	= P->p["-o"]+ job_name+ ".prelim_bidir_hits.bed";
		model_run(P, rank, nprocs,0, job_ID, LG);
		
	}
	LG->write("exiting bidir module....................................done\n\n", verbose);
	return 1;
}
