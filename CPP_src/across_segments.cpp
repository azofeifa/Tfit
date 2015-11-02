#include <mpi.h>
#include "load.h"
#include "across_segments.h"
#include "model.h"
#include "template_matching.h"
#include "model_single.h"
#include <iostream>
#include <fstream>
#include <map>
#include <time.h>
#include "omp.h"
#include "read_in_parameters.h"
using namespace std;


string check_file(string FILE, int i){ //don't want to write over an existing file
	string template_file 	= FILE + to_string(i);
	ifstream FH(template_file);
	if (FH){
		FH.close();
		return	check_file(FILE, i+1);
	}
	return template_file;
}

map<int,vector<classifier> > initialize_data_struct(int minK, int maxK, 
	int rounds, int num_proc, double scale, 
	double max_noise,  double convergence_tresh, 
	int max_iterations, double r_mu, double ALPHA_0, double BETA_0, 
	double ALPHA_1, double BETA_1, double ALPHA_2, double ALPHA_3){
	map<int,vector<classifier> >  data_struct;	
	for (int k = minK; k<=maxK;k++){
		data_struct[k] 	= vector<classifier>(rounds);
		for (int r = 0; r<rounds; r++){
			data_struct[k][r] 	= classifier(k,  convergence_tresh, max_iterations, 
					max_noise,r_mu, ALPHA_0, BETA_0, ALPHA_1, BETA_1, ALPHA_2, ALPHA_3,0);
		}
	}
	return data_struct;
}

map<int, classifier> getMax(map<int,vector<classifier> > DS){
	typedef map<int, vector<classifier>>::iterator it_type;
	map<int, classifier> reduced;
	double MAX=nINF;
	for(it_type K = DS.begin(); K != DS.end(); K++) {
		for (int i = 0; i < K->second.size(); i ++){
			if (DS[K->first][i].ll > MAX and DS[K->first][i].converged) {
				reduced[K->first] 	= DS[K->first][i];
				MAX 				= DS[K->first][i].ll;
			}
		}
	}
	return reduced;

}


string get_header(params * P){
	string header 	= currentDateTime() + ",";
	header+="minK: " + P->p["-minK"] + ",";
	header+="maxK: " +  P->p["-maxK"] + ",";
	header+="scale: " +  P->p["-ns"] + ",";
	header+="bin_resolution: " + P->p["-br"] + ",";
	header+="convergence_tresh: " + P->p["-ct"] + ",";
	header+="max_iterations: " + P->p["-mi"] + ",";
	header+="spec_chrom: " + P->p["-chr"] + ",";
	header+="max_noise: " + P->p["-max_noise"] + ",";
	header+="r_mu: " + P->p["-r_mu"] +",";
	header+="rounds: " + P->p["-rounds"]+",";
	header+="ALPHA_0: " + P->p["-ALPHA_0"] +",";
	header+="BETA_0: " + P->p["-BETA_0"] +",";
	header+="ALPHA_1: " + P->p["-ALPHA_1"] +",";
	header+="BETA_1: " + P->p["-BETA_1"] +",";
	header+="ALPHA_2: " + P->p["-ALPHA_2"] +",";
	header+="ALPHA_3: " + P->p["-ALPHA_3"]  ;
	
	header+="\n";
	return header;
}

void run_model_accross_segments(vector<segment*> segments, params *P){
	typedef map<int, classifier>::iterator it;
	



	int N 			= segments.size();
	string out_dir 	= P->p["-o"]; 
	string out_file_template 	= out_dir+"model_fits_out_" + P->p["-chr"]+"_";
	string out_file 			= check_file(out_file_template, 1);
	
	

	ofstream FHW;
	
	FHW.open(out_file);

	
	FHW<<get_header(P);

	int minK 					= stoi(P->p["-minK"]);
	int maxK 					= stoi(P->p["-maxK"]);
	double max_noise 			= stod(P->p["-max_noise"]);
	double convergence_tresh 	= stod(P->p["-ct"]);
	int max_iterations 			= stoi(P->p["-mi"]);
	double r_mu 				= stod(P->p["-r_mu"]);
	int num_proc 				= stoi(P->p["-np"]);
	int rounds 					= stoi(P->p["-rounds"]);
	double scale 				= stod(P->p["-ns"]);
	double ALPHA_0 				= stod(P->p["-ALPHA_0"]);
	double BETA_0 				= stod(P->p["-BETA_0"]);
	double ALPHA_1 				= stod(P->p["-ALPHA_1"]);
	double BETA_1 				= stod(P->p["-BETA_1"]);
	double ALPHA_2 				= stod(P->p["-ALPHA_2"]);
	double ALPHA_3 				= stod(P->p["-ALPHA_3"]);
	
	for (int i = 0; i < N; i++ ){
		if (segments[i]->N > 0){
			vector<double> mu_seeds 		=  peak_bidirs(segments[i]);
			map<int,vector<classifier> > DS = initialize_data_struct(minK, maxK, 
				rounds, num_proc, scale, max_noise,  
				convergence_tresh, max_iterations,r_mu, ALPHA_0,  BETA_0, 
	 		ALPHA_1,  BETA_1,  ALPHA_2,  ALPHA_3);
			

			classifier clf(0, convergence_tresh, max_iterations, 
						max_noise,r_mu, ALPHA_0, BETA_0, ALPHA_1, BETA_1, ALPHA_2, ALPHA_3,0);
			clf.fit(segments[i], mu_seeds);
			
			FHW<<segments[i]->write_out();
			FHW<<"~0"<<","<<to_string(clf.ll)<<",1,0"<<endl;
			FHW<<"U: "<<to_string(segments[i]->minX)<<","<<to_string(segments[i]->maxX)<<",1,"<<to_string(clf.pi)<<endl;
			for (int k = minK; k <=maxK;k++ ){
				#pragma omp parallel for num_threads(num_proc)
				for (int j = 0; j < rounds; j++){
					DS[k][j].fit(segments[i], mu_seeds);
				}
			}
			//*********************************************
			//we want to output for each model, k, the highest log likelihood estimate
			//and then the respective parameter estimates for each component
			map<int, classifier> reduced 	= getMax(DS);
			//okay now lets write it out to the file
			for (int k = minK; k<=maxK; k++){
				if (reduced.find(k) !=reduced.end()){
					classifier K 	= reduced[k];
					FHW<<"~"<<to_string(k)<<","<<to_string(K.ll)<<","<<to_string(K.converged)<<","<<to_string(K.last_diff)<<endl;
					FHW<<K.print_out_components();
				}else{
					FHW<<"~"<<to_string(k)<<","<<to_string(nINF)<<","<<to_string(0)<<","<<to_string(INF)<<endl;
				}
			}
		}		
	}
}



void free_segments(vector<segment*> segments){
	for (int i = 0; i < segments.size(); i++){
		delete segments[i];
	}
}



vector<classifier> get_vector_classifiers(params * P, int K){
	vector<classifier> clfs(stoi(P->p["-rounds"]));
	for (int i =0 ; i< stoi(P->p["-rounds"]); i++ ){
		clfs[i] 	= classifier(K, stod(P->p["-ct"]), stoi(P->p["-mi"]), stod(P->p["-max_noise"]), 
			stod(P->p["-r_mu"]), stod(P->p["-ALPHA_0"]), stod(P->p["-BETA_0"]), stod(P->p["-ALPHA_1"]), 
			stod(P->p["-BETA_1"]), stod(P->p["-ALPHA_2"]) , stod(P->p["-ALPHA_3"]),0 );
	}
	return clfs;
}

vector<classifier> get_vector_classifiers2(params * P, int K){

	int res 	= stoi(P->p4["-foot_res"]);
	double lower, upper;
	lower=0, upper=1000;

	double delta; 
	if (res==0){
		delta 		= 0;
	}else{
		delta 	= (upper-lower) / float(res);
	}
	vector<classifier> clfs(stoi(P->p4["-rounds"])*res   );
	double foot_print;
	int i 	= 0;
	int r 	= 0;
	double scale 	= stod(P->p4["-ns"]);
	while (i < clfs.size()){
		if (r > res){
			r=0;
		}
		foot_print 	= (delta*r + lower)/scale;
		clfs[i] 	= classifier(K, stod(P->p4["-ct"]), stoi(P->p4["-mi"]), stod(P->p4["-max_noise"]), 
			stod(P->p4["-r_mu"]), stod(P->p4["-ALPHA_0"]), stod(P->p4["-BETA_0"]), stod(P->p4["-ALPHA_1"]), 
			stod(P->p4["-BETA_1"]), stod(P->p4["-ALPHA_2"]) , stod(P->p4["-ALPHA_3"]), false,foot_print );
		i++;
		r++;
	
	}
	return clfs;
}
classifier fit_noise(params *P){
	classifier clf(0, stod(P->p["-ct"]),
					stoi(P->p["-mi"]), 
					stod(P->p["-max_noise"]), 
					stod(P->p["-r_mu"]), 
					stod(P->p["-ALPHA_0"]), 
					stod(P->p["-BETA_0"]), 
					stod(P->p["-ALPHA_1"]), 
					stod(P->p["-BETA_1"]), 
					stod(P->p["-ALPHA_2"]), 
					stod(P->p["-ALPHA_3"]),0 									
					);
	return clf;
}
vector<simple_c> get_max(vector<classifier> clfs, 
	double noise_ll, int seg, int complexity, int found_bidirs,
	 int bidir_ID, double NN){
	double max = nINF;
	classifier argmax;
	vector<simple_c> scs;
	simple_c sc;
	for (int i = 0; i < clfs.size(); i++){
		if (clfs[i].ll > max and clfs[i].converged){
			max=clfs[i].ll, argmax=clfs[i];	
		}
	}
	if (max==nINF){	
		sc.IDS[0] 	= seg, sc.IDS[1] 	= complexity;
		sc.IDS[2] 	= found_bidirs, sc.IDS[3] = bidir_ID;
		sc.ps[0] 	= 0;
		sc.ps[1] 	= 0;
		sc.ps[2] 	= 0;
		sc.ps[3] 	= 0;
		sc.ps[4] 	= 0;
		sc.ps[5] 	= 0; 
		sc.ps[6] 	= 0;
		sc.ps[7] 	= 0; 
		sc.ps[8] 	= 0;
		sc.ps[9] 	= 0;
		sc.ps[10] 	= 0;
		sc.ps[11] 	= NN;
		sc.ps[12] 	= 0;
		sc.ll 		= max, sc.noise_ll 	= noise_ll;		
		scs.push_back(sc);
	}else{
		for (int c = 0; c < argmax.K; c++){
			simple_c scc;
	
			scc.IDS[0] 	= seg, scc.IDS[1] 	= complexity;
			scc.IDS[2] 	= found_bidirs, scc.IDS[3] = bidir_ID;
			scc.ps[0] 	= argmax.components[c].bidir.mu;
			scc.ps[1] 	= argmax.components[c].bidir.si;
			scc.ps[2] 	= argmax.components[c].bidir.l;
			scc.ps[3] 	= argmax.components[c].bidir.w;
			scc.ps[4] 	= argmax.components[c].bidir.pi;
			scc.ps[5] 	= argmax.components[c].forward.w; 
			scc.ps[6] 	= argmax.components[c].reverse.w;
			scc.ps[7] 	= argmax.components[c].forward.b; 
			scc.ps[8] 	= argmax.components[c].reverse.a;
			scc.ps[9] 	= argmax.components[c].forward.pi; 
			scc.ps[10] 	= argmax.components[c].reverse.pi;
			
			scc.ps[11] 	= NN;
			scc.ps[12] 	= argmax.components[c].bidir.foot_print;
			scc.ll 	= max, scc.noise_ll 	= noise_ll;		

			scs.push_back(scc);
			
		}
		
	}
	return 	scs;
	
}





vector<simple_c> wrapper_pp(segment * s, params * P, int seg){
	
	int num_proc 				= stoi(P->p["-np"]);
	double scale 				= stod(P->p["-ns"]);
	
	s->insert_bidirectional_data(
			int(stoi(P->p["-pad"])/stod(P->p["-ns"]) ));
	vector<simple_c> fits;
	for (int j =0; j < s->bidirectional_data.size(); j++){
		s->bidirectional_data[j]->bin(1  , 0, true);
		s->bidirectional_data[j]->SCALE 	= scale;
		vector<double> mu_seeds 			= peak_bidirs(s->bidirectional_data[j]);
		classifier noise_clf(0, stod(P->p["-ct"]), stoi(P->p["-mi"]), stod(P->p["-max_noise"]), 
			stod(P->p["-r_mu"]), stod(P->p["-ALPHA_0"]), stod(P->p["-BETA_0"]), stod(P->p["-ALPHA_1"]), 
			stod(P->p["-BETA_1"]), stod(P->p["-ALPHA_2"]) , stod(P->p["-ALPHA_3"]),0 );
		noise_clf.fit(s->bidirectional_data[j], mu_seeds);
		double noise_ll 	= noise_clf.ll;
		for (int k =1 ; k <= s->bidirectional_data[j]->counts; k++ ){
			vector<classifier> 	clfs 			= get_vector_classifiers(P, 
				k);
			#pragma omp parallel for num_threads(num_proc)
			for (int t = 0; t <  stoi(P->p["-rounds"]); t++){
				clfs[t].fit(s->bidirectional_data[j], mu_seeds);
			}
			vector<simple_c> bidir_components = get_max(clfs, noise_ll, seg, k, 
				s->bidirectional_data[j]->counts, j, s->bidirectional_data[j]->N );
			for (int b=0; b < bidir_components.size(); b++){
				fits.push_back(bidir_components[b]);
			}
			
		}
	}
	return fits;
}

vector<simple_c> wrapper_pp_just_segments(segment * s , params * P, int seg, int np){
	vector<simple_c> fits;
	classifier noise_clf(0, stod(P->p4["-ct"]), stoi(P->p4["-mi"]), stod(P->p4["-max_noise"]), 
		stod(P->p4["-r_mu"]), stod(P->p4["-ALPHA_0"]), stod(P->p4["-BETA_0"]), stod(P->p4["-ALPHA_1"]), 
		stod(P->p4["-BETA_1"]), stod(P->p4["-ALPHA_2"]) , stod(P->p4["-ALPHA_3"]), 0. );
	
	noise_clf.fit(s, s->centers);
	
	double noise_ll 	= noise_clf.ll;
	for (int k = s->counts; k<= s->counts+1;k++){
		vector<classifier> 	clfs 			= get_vector_classifiers2(P,k);
		#pragma omp parallel for num_threads(np)
		for (int t = 0; t <  clfs.size(); t++){
			clfs[t].fit(s, s->centers);
		}
		vector<simple_c> bidir_components = get_max(clfs, noise_ll, seg, k, s->counts, 1, s->N );
		for (int b=0; b < bidir_components.size(); b++){
			fits.push_back(bidir_components[b]);
		}
		
	}
	
	return fits;
}

vector<simple_c> bidir_components_to_simplec(vector<classifier> clfs, vector<segment *> FSI){
	vector<simple_c> fits;
	typedef vector<classifier>::iterator it_type;
	int i 	= 0;
	for (it_type c=clfs.begin(); c!=clfs.end(); c++){
		double ll 	= (*c).ll;
		double NN 	= FSI[i]->N;
		for (int k = 0 ; k < (*c).init_parameters.size(); k++ ){
			simple_c scc;
			scc.IDS[0] 	= i, scc.IDS[1] 	= (*c).init_parameters.size();
			scc.IDS[2] 	= FSI[i]->ID, scc.IDS[3] 	= k;
			scc.ps[0] 	= (*c).components[k].bidir.mu;
			scc.ps[1] 	= (*c).components[k].bidir.si;
			scc.ps[2] 	= (*c).components[k].bidir.l;
			scc.ps[3] 	= (*c).components[k].bidir.w;
			scc.ps[4] 	= (*c).components[k].bidir.pi;
			scc.ps[5] 	= (*c).components[k].forward.w; 
			scc.ps[6] 	= (*c).components[k].reverse.w;
			scc.ps[7] 	= (*c).components[k].forward.b; 
			scc.ps[8] 	= (*c).components[k].reverse.a;
			scc.ps[9] 	= (*c).components[k].forward.pi; 
			scc.ps[10] 	= (*c).components[k].reverse.pi;
			scc.ps[11] 	= NN;
			scc.ll 		= ll, scc.noise_ll 	= nINF;		
			fits.push_back(scc);
		}
		i++;
	}

	return fits;
}


vector<simple_c> run_model_accross_segments_template(vector<segment*> segments, 
	params *P){
	vector<simple_c> fits;
	double noise_LL, LL, mu,si, l, w_e, w_f, w_r, pi, pi_f, pi_r, r_a, f_b;
	for (int i = 0; i < segments.size(); i++){
		vector<simple_c> curr_fits 		= wrapper_pp(segments[i], P, i); 
    	for (int j = 0; j < curr_fits.size(); j++){
    		fits.push_back(curr_fits[j]);
    	}
    }
    //we want to return a vector of very very simple structs
    return fits;
}

vector<simple_c> run_model_accross_segments_to_simple_c(vector<segment *> segments, 
	params * P, ofstream& log_file){
	vector<simple_c> fits;
	log_file<<"(across_segments) model progress...";
	log_file<<to_string(0) + "%%,";
	log_file.flush();

	double percent 	= 0;
	double N 		= segments.size();
	int rounds 		= stod(P->p4["-rounds"]);
	int inner_np 	= 1;
	inner_np 	= omp_get_max_threads();
	
	map<int, vector<simple_c>> 	G;
	int PP 	= 0;
	for (int i = 0; i < segments.size(); i++){
		if ((i/N) > percent+ 0.1 ){
			PP+=10;
			log_file<<to_string(PP) + "%%,";
			log_file.flush();
			percent 	= i / N;
		}
		vector<simple_c> curr_fits 	= wrapper_pp_just_segments(segments[i], P , i, inner_np);
		G[i] 	= curr_fits;
	}
	typedef map<int, vector<simple_c> >::iterator it_type;
	for (it_type i 	= G.begin(); i!=G.end(); i++){
		for (int j = 0; j < i->second.size(); j++){
			fits.push_back(i->second[j]);
	    }
    }
	log_file<<"...done\n";
	log_file.flush();
	return fits;
}

vector<simple_c> move_elongation_support(vector<segment *> FSI, params * P){
	int num_proc 				= omp_get_max_threads();
	vector<classifier> clfs(int(FSI.size() ));
	vector<simple_c> fits;
	int rounds 					= 5;
	#pragma omp parallel for num_threads(num_proc)
	for (int i = 0; i < FSI.size(); i++){
		classifier clf(stod(P->p4["-ct"]), stoi(P->p4["-mi"]), stod(P->p4["-max_noise"]), 
			stod(P->p4["-r_mu"]), stod(P->p4["-ALPHA_0"]), stod(P->p4["-BETA_0"]), stod(P->p4["-ALPHA_1"]), 
			stod(P->p4["-BETA_1"]), stod(P->p4["-ALPHA_2"]) , stod(P->p4["-ALPHA_3"]) ,FSI[i]->fitted_bidirs ,0 ) 	 ;
		clf.fit_uniform_only2(FSI[i]);

	
	}
	//convert to simple_c
	fits 	= bidir_components_to_simplec(clfs, FSI);
	return fits;
}


double BIC_score(double ll, double N, int type){
	if (type == 0){
		return -2*ll + 1*log(N);
	}else if (type==3){
		return -2*ll + 6*log(N);		
	}else{
		return -2*ll + 3*log(N);			
	}
}
single_simple_c classifier_single_to_simple_c(NLR arg_clf, double maxll, segment *data, int K, double NN){
	single_simple_c sc; 
	for (int i = 0; i< 5; i++){
		if (i < data->chrom.size()){
			sc.chrom[i] 	= data->chrom[i];
		}else{
			sc.chrom[i] 	= '\0';
		}
	}
	sc.st_sp[0] 	= data->start, sc.st_sp[1]=data->stop, sc.st_sp[2]=data->ID;
	sc.ps[0] 	= arg_clf.mu, sc.ps[1] 	= arg_clf.si, sc.ps[2] = arg_clf.wn; 
	sc.ps[3] 	= arg_clf.l, sc.ps[4] = arg_clf.r, sc.ps[5] = arg_clf.wl, sc.ps[6] = arg_clf.wr;
	sc.ps[7] 	= maxll, sc.ps[8] 	= NN;
	return sc;

}

single_simple_c classifier_single_to_simple_c_noise(segment* data, int i, double ll, double NN){

	single_simple_c sc; 
	for (int i = 0; i< 5; i++){
		if (i < data->chrom.size()){
			sc.chrom[i] 	= data->chrom[i];
		}else{
			sc.chrom[i] 	= '\0';
		}
	}
	sc.st_sp[0] 	= data->start, sc.st_sp[1]=data->stop, sc.st_sp[2]=data->ID;
	sc.ps[0] 	= 0, sc.ps[1] 	= 0, sc.ps[2] = 0; 
	sc.ps[3] 	= data->minX, sc.ps[4] = data->maxX, sc.ps[5] = 1.0, sc.ps[6] = 1.0;
	sc.ps[7] 	= ll, sc.ps[8] = NN;
	return sc;
			
}


vector<single_simple_c> run_single_model_across_segments(vector<segment *> FSI, params * P){
	int N 			= FSI.size();
	int num_proc 	= stoi(P->p5["-np"]);
	int rounds 		= stoi(P->p5["-rounds"]);
	double scale 	= stod(P->p5["-ns"]);
	int mi 	= stod(P->p5["-mi"]);
	double ct 	= stod(P->p5["-ct"]);
	vector<single_simple_c> fits;
	vector<vector<classifier_single>> clf_fits(N);
	#pragma omp parallel for num_threads(num_proc)	
	for (int i 	= 0; i < N;i++){
		//want to fit four models, noise, left, right, both
		double BIC_min 	= INF;

		int arg_type;
		double ll, best_ll;
		classifier_single BIC_best;
		vector<classifier_single> current;
		for (int j = 0; j < 4; j++){
			double maxll 	= nINF;
			classifier_single  arg_clf;
			if (j > 0){
				bool SET 	= false;
				for (int r = 0; r < rounds; r++){
					classifier_single clf(ct, mi, 1,j, scale);
					ll 				= clf.fit(FSI[i]);
					if (ll > nINF){
						maxll 		= ll;
						arg_clf 	= clf;
						SET 		= true;
					}else if(r == rounds-1 and not SET){
						maxll 		= ll;
						arg_clf 	= clf;	
					}
				}
			}else{
				classifier_single clf(ct, mi, 0,j, scale);	
				maxll 		= clf.fit(FSI[i]);
				arg_clf 	= clf;
			}
			//calc BIC score
			current.push_back(arg_clf);
			if (BIC_score(maxll, FSI[i]->N, j) <  BIC_min  ){

				BIC_min 	= BIC_score(maxll, FSI[i]->N, j), arg_type = j;
				BIC_best 	= arg_clf;		
				best_ll 	= maxll;
			}

		}
		//transform to simple_c;
		clf_fits[i] = current;
		
		
	}

	for (int i = 0 ; i < N; i++){
		for (int j = 0 ; j < clf_fits[i].size(); j++)
			if (clf_fits[i][j].K==0){
				fits.push_back(classifier_single_to_simple_c_noise(FSI[i], i, clf_fits[i][j].ll,FSI[i]->N ));
			}else{
				for (int c = 0; c < clf_fits[i][j].K; c++){
					fits.push_back(classifier_single_to_simple_c(clf_fits[i][j].components[c],clf_fits[i][j].ll, FSI[i],clf_fits[i][j].K ,FSI[i]->N ));
				}
			}

	}
	return fits;
}

map<int, vector<classifier> > make_classifier_struct_free_model(params * P, segment * data){

	int min_k 	= stoi(P->p["-minK"]);
	int max_k 	= stoi(P->p["-maxK"]);
	int rounds 	= stoi(P->p["-rounds"]);
	int BDS 	= int(data->centers.size());
	map<int, vector<classifier> > A;
	for (int k =max(min_k, BDS/2); k < BDS; k++){
		for (int r = 0; r < rounds; r++){
			A[k].push_back(classifier(k, stod(P->p["-ct"]), stoi(P->p["-mi"]), stod(P->p["-max_noise"]), 
			stod(P->p["-r_mu"]), stod(P->p["-ALPHA_0"]), stod(P->p["-BETA_0"]), stod(P->p["-ALPHA_1"]), 
			stod(P->p["-BETA_1"]), stod(P->p["-ALPHA_2"]) , stod(P->p["-ALPHA_3"]),0 ));
		}
	}
	return A;

}


simple_c_free_mode::simple_c_free_mode(bool FOUND, double ll, 
	component C, int K, segment * data, int i, double forward_N, double reverse_N){
	SS[0]=ll, SS[1]=forward_N, SS[2]=reverse_N;
	ID[0]=data->ID, ID[1]=data->start, ID[2]=data->stop, ID[3]=K;
	if (FOUND){
		ID[4]=1; 
	}else{
		ID[4]=0;//didnt converge
	}
	for (int c =0 ; c < 5; c++){
		if (c < data->chrom.size()){
			chrom[c]=data->chrom[c];
		}else{
			chrom[c]='\0';
		}
	}
	chrom[5]='\0';
	if (not FOUND){
		for (int i = 0; i < 11; i++ ){
			ps[i] 	= 0;
		}
	}else{
		ps[0]=C.bidir.mu,ps[1]=C.bidir.si,ps[2]=C.bidir.l, ps[3]=C.bidir.w, ps[4]=C.bidir.pi;
		ps[5]=C.forward.b, ps[6]=C.forward.w, ps[7]=C.forward.pi;
		ps[8]=C.reverse.a, ps[9]=C.reverse.w, ps[10]=C.reverse.pi;
	}
}
simple_c_free_mode::simple_c_free_mode(){}

vector<simple_c_free_mode> transform_free_mode(bool FOUND, double ll, component * components, 
	int K, segment * data, int i, double forward_N, double reverse_N) {
	vector<simple_c_free_mode> SC(K);
	
	for (int k = 0 ; k < K; k++){
		SC[k]=simple_c_free_mode(FOUND, ll, components[k], K, data, i,forward_N, reverse_N);
	}
	return SC;

	
}

map<int, vector<simple_c_free_mode> > get_max_from_free_mode(map<int, vector<classifier> > A, segment * data, int i){
	map<int, vector<simple_c_free_mode> > BEST;
	typedef map<int, vector<classifier> >::iterator it_type_A;
	//get forward and reverse N
	double forward_N =0, reverse_N=0;
	for (int i  = 0 ; i < data->XN;i++){
		forward_N+=data->X[1][i];
		reverse_N+=data->X[2][i];
	}


	for (it_type_A a = A.begin(); a!=A.end(); a++){
		component * best_components;
		double best_ll 	= nINF;
		bool FOUND 		= false;
		for (int r = 0; r < a->second.size(); r++){
			if (A[a->first][r].ll > best_ll){
				best_ll 			= A[a->first][r].ll;
				best_components 	= A[a->first][r].components;
				FOUND 				= true;
			}
		}
		BEST[a->first] 	= transform_free_mode(FOUND, best_ll, best_components, a->first, data, i, forward_N, reverse_N);

	}

	return BEST;
}


vector<map<int, vector<simple_c_free_mode> >> run_model_across_free_mode(vector<segment *> FSI, params * P, ofstream& log_file){
	vector<map<int, vector<simple_c_free_mode> >> D;
	typedef map<int, vector<classifier> > ::iterator it_type;
	double scale 	= stof(P->p["-ns"]);
	int num_proc 				= omp_get_max_threads();
	log_file<<"(across_segments) running model across provided intervals...";
	log_file.flush();
	double N 		= FSI.size();
	double percent 	= 0;
	int elon_move 	= stoi(P->p["-elon"]);
	int topology 	= stoi(P->p["-topology"]);
	for (int i = 0 ; i < FSI.size(); i++){
		if ((i / N) > (percent+0.05)){
			log_file<<to_string(int((i / N)*100))+"%,";
			log_file.flush();
			percent 	= (i / N);
		}

		//first need to populate data->centers
		for (int b = 0 ; b < FSI[i]->bidirectional_bounds.size(); b++){
			double center = FSI[i]->bidirectional_bounds[b][0] +  FSI[i]->bidirectional_bounds[b][1] ;
			center/=2.;
			center-=FSI[i]->start;
			center/=scale;
			 FSI[i]->centers.push_back(center);
		}
		segment * data 	= FSI[i];
		map<int, vector<classifier> > A 	= make_classifier_struct_free_model(P, FSI[i]);
		for (it_type k = A.begin(); k!= A.end(); k++){
			int N 	=  k->second.size();
			#pragma omp parallel for num_threads(num_proc)	
			for (int r = 0; r < N; r++ ){
				A[k->first][r].fit2(data, data->centers,topology,elon_move);
			}
		}
		D.push_back(get_max_from_free_mode(A, FSI[i], i));
	}
	log_file<<"100%";
	log_file.flush();
	return D;
}























