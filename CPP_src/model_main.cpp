#include "model_main.h"
#include "template_matching.h"
#include "density_profiler.h"
#include "MPI_comm.h"
#include <omp.h>
using namespace std;
int model_run(params * P, int rank, int nprocs, double density, int job_ID, Log_File * LG){
	int verbose 	= stoi(P->p["-v"]);
	LG->write("\ninitializing model module...............................done\n\n",verbose);
	int threads 	= omp_get_max_threads();//number of openMP threads that are available for use	
	string job_name = P->p["-N"];

	//=======================================================================================
	//input file paths
	string forward_bed_graph_file 	= P->p["-i"];
	string reverse_bed_graph_file 	= P->p["-j"];
	string interval_file 			= P->p["-k"];
	string out_file_dir 			= P->p["-o"];
	string spec_chrom 				= P->p["-chr"];
	//=======================================================================================
	//(1a) load intervals and keep track of their associated IDS
	map<int, string> IDS;
	vector<segment *> FSI;
	LG->write("loading intervals of interest...........................",verbose);
	if (rank==0){
		FSI 	= load::load_intervals_of_interest(interval_file, IDS, P );
	}
	LG->write("done\n",verbose);
	//(1b) now broad cast out the intervals of interest to individual MPI processes
	LG->write("sending interval assignments............................",verbose);
	map<string, vector<segment *> > GG 	= MPI_comm::send_out_single_fit_assignments(FSI, rank, nprocs);
	LG->write("done\n",verbose);


	//=======================================================================================
	//(2a) load bedgraph files and insert them into intervals of interest (interval tree...)
	LG->write("inserting bedgraph data.................................",verbose);
	vector<segment*> integrated_segments= load::insert_bedgraph_to_segment_joint(GG, 
		forward_bed_graph_file, reverse_bed_graph_file, rank);
	//(2b) for each segment we are going to bin and scale and center, numerical stability
	LG->write("done\n",verbose);
	LG->write("binning, centering, scaling.............................",verbose);
	load::BIN(integrated_segments, stod(P->p["-br"]), stod(P->p["-ns"]),true);	
	LG->write("done\n",verbose);

	//=======================================================================================
	//(3a) now run template matching for seeding the EM  
	LG->write("running template matching...............................",verbose);
	run_global_template_matching(integrated_segments, out_file_dir, 3, 
			0.0,stod(P->p["-ns"]),stod(P->p["-bct"]), threads,0. ,0);	
	LG->write("done\n",verbose);
	//=======================================================================================
	//(4a) now going to run the model across all segments

	vector<map<int, vector<simple_c_free_mode> >> FITS 		= run_model_across_free_mode(integrated_segments,
		 P,LG);
	//(4b) gather all the model fits
	LG->write("gathering all model fits................................",verbose);
	map<int, map<int, vector<simple_c_free_mode>  > > GGG 	= MPI_comm::gather_all_simple_c_free_mode(FITS, rank, 
		nprocs);
	LG->write("done\n",verbose);
	//(4c) now write out model fits
	if (rank==0){//write_out_to_MLE, //out_dir+  P->p["-N"] + "-" + to_string(job_ID)+  "_K_models_MLE.tsv"
		LG->write("writing out results (MLE)...............................",verbose);
		string file_name = "";
		load::write_out_models_from_free_mode(GGG, P, job_ID, IDS, density, file_name);
		LG->write("done\n",verbose);
		LG->write("loading results (MLE)...................................",verbose);
		vector<segment_fits *> fits 		= load::load_K_models_out(file_name);
		LG->write("done\n",verbose);		
		LG->write("writing out results (model selection)...................",verbose);
		load::write_out_bidirectionals_ms_pen(fits, P, job_ID, density );
		LG->write("done\n",verbose);
	}
	LG->write("\nexiting model module....................................done\n\n",verbose);
	//wait for everybody to catch up
	MPI_comm::wait_on_root(rank, nprocs);



	return 1;




}
