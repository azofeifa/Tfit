#include "bootstrap.h"
#include "load.h"
#include <vector>
#include "read_in_parameters.h"
#include "model.h"
#include <random>
#include <iostream>
#include <fstream>

#include <omp.h>
#include "template_matching.h"
using namespace std;

int sample(double ** CDF, int XN, double sum_N, segment * NS, double pi , segment * S){
	random_device rd;
	mt19937 MT(rd());
	
	uniform_real_distribution<double> RAND(0,1);
	double  U;  
	int i 	= 0;
	int ct 	= 0;
	int j 	= 0;
	int s 	= 0;
	while (ct < sum_N){
		U 		= RAND(MT);
		if (U < pi ){
			s 	= 1;
		}else{
			s 	= 2;
		}
		j 		= 0;
		U 		= RAND(MT);
		while (j+1 < XN and CDF[s][j] < U){
			j++;
		}
		if (S->X[s][j] > 0){
			NS->X[s][j]+=1.;
			ct++;	
		}

	}
	

	return i;
	
}

void subsample(segment * S, segment * NS ){
	double ** CDF 	= new double*[3];
	NS->X 			= new double*[3];
	NS->minX = S->minX, NS->maxX = S->maxX;
	NS->XN 			= S->XN;
	NS->N 			= S->N;
	NS->SCALE 		= S->SCALE;
	int BINS 		= int(S->XN);
	for (int j = 0; j < 3; j++){
		CDF[j]=new double[BINS], NS->X[j]=new double[BINS];
	}
	for (int i = 0 ; i< S->XN; i++){
		CDF[0][i] 	= S->X[0][i], NS->X[0][i] = S->X[0][i];
		CDF[1][i] 	= 0,CDF[2][i] 	= 0;
		NS->X[1][i] = 0,NS->X[2][i] = 0;		
	}
	double forward_sum =0, reverse_sum = 0, sum_N = 0, pi = 0;
	for (int i = 0; i < S->XN; i++){
		forward_sum+=S->X[1][i];
		reverse_sum+=S->X[2][i];
	}
	sum_N 	= forward_sum+reverse_sum;
	double forward_N = forward_sum;
	double reverse_N = reverse_sum;
	pi 		= forward_sum / sum_N;
	forward_sum = 0, reverse_sum = 0;
	for (int i = 0; i < S->XN; i++){
		forward_sum+=S->X[1][i];
		reverse_sum+=S->X[2][i];
		CDF[1][i ] 	= forward_sum / forward_N;
		CDF[2][i ] 	= reverse_sum / reverse_N;
	}

	sample( CDF, S->XN, S->N, NS , pi, S);
	
}

vector<segment *> make_bootstraps(segment * S, params * P){
	int rounds 	= stoi(P->p6["-rounds"]);
	int brounds = stoi(P->p6["-brounds"]);
	int np 		= omp_get_max_threads();
	//make segmnets
	vector<segment *> segments(brounds);
	// NS->minX = S->minX, NS->maxX = S->maxX;
	// NS->XN 			= S->XN;
	// NS->N 			= S->N;
	// NS->SCALE 		= S->SCALE;
	
	#pragma omp parallel for num_threads(np)
	for (int r = 0; r < brounds; r++){
		//subsample get new segment
		segment * T 	= S;
		segment * NS 	= new segment(S->chrom, S->start, S->stop);
		NS->parameters 	= S->parameters;

		subsample(  S,   NS );
		segments[r] 	= NS;
	}
	return segments;
}


vector<vector<double>> sort_bootstrap_parameters(vector<vector<double>> X){
	bool changed=true;
	while (changed){
		changed=false;
		for (int i = 0; i < X.size()-1; i++  )	{
			if (X[i][0] > X[i+1][0]){ //sort by starting position
				vector<double> copy 	= X[i];
				X[i] 					= X[i+1];
				X[i+1] 					= copy;
				changed=true;
			}
		}
	}
	return X;	
}

double get_mean(vector<double> X){
	double S 	= 0;
	double N 	= 0;
	for (int i = 0; i < X.size(); i++){
		S+=X[i];
		N+=1;
	}
	return S/N;
}

void run_bootstrap_across(vector<segment *> segments, params * P, ofstream& log_file){
	int rounds 		= stoi(P->p6["-rounds"]);
	double scale 	= stod(P->p6["-ns"]);
	int np 			= omp_get_max_threads();
	int res 	= stoi(P->p6["-foot_res"]);
	double lower=0, upper=500, delta= 0;
	log_file<<"(across_segments) model progress...";
	log_file<<to_string(0) + "%%,";
	log_file.flush();

	if (res!=0){
		delta 	= (upper-lower) / float(res);
	}
	double NN 		= segments.size();
	double percent  = 0.;
	int PP 			= 0;
	for (int i = 0; i < segments.size(); i++){
		if ((i/NN) > percent+ 0.1 ){
			PP+=10;
			log_file<<to_string(PP) + "%%,";
			log_file.flush();
			percent 	= i / NN;
		}
		
		int K 				= segments[i]->parameters.size();
		segments[i]->parameters 	= sort_bootstrap_parameters(segments[i]->parameters);
		vector<double> mu_seeds( K ); 

		for (int c = 0; c < segments[i]->parameters.size(); c++){
			mu_seeds[c]=(segments[i]->parameters[c][0] - segments[i]->start) / scale;
		}
		map<int, map<int, vector<double> >> variances;
		map<int,  vector<double> >  final_variances;
		if (segments[i]->N > 0){
			vector<segment *> bootstraps = make_bootstraps(segments[i], P);
			int B 				= bootstraps.size();
			vector<classifier> clfs( B );
			vector<classifier> fits(B);
			double foot_print;
			#pragma omp parallel for num_threads(np)
			for (int b = 0; b < bootstraps.size(); b++ ){
				double best_ll 	= nINF;
				classifier best_clf;
				segment * NS 	= bootstraps[b];
				for (int w = 0; w < res;w++){
					foot_print 	= (delta*w + lower)/scale;
					for (int r = 0 ; r < rounds; r++ ){
						classifier current_clf(K, stod(P->p6["-ct"]), stoi(P->p6["-mi"]), stod(P->p6["-max_noise"]), 
							stod(P->p6["-r_mu"]), stod(P->p6["-ALPHA_0"]), stod(P->p6["-BETA_0"]), stod(P->p6["-ALPHA_1"]), 
							stod(P->p6["-BETA_1"]), stod(P->p6["-ALPHA_2"]) , stod(P->p6["-ALPHA_3"]), false,foot_print );
						current_clf.fit(NS, mu_seeds);
						if (w == 0 or current_clf.ll > best_ll){
							best_clf 	= current_clf;
							best_ll 	= current_clf.ll;
						}
					}
				}
				fits[b]=best_clf;
			}
			for (int b =0 ; b < B; b++){
				vector<vector<double>> bootstrapped_parameters(K);
				for (int k = 0; k < K; k++){
					vector<double> current_parameters(6);
					current_parameters[0] = fits[b].components[k].bidir.mu;
					current_parameters[1] = fits[b].components[k].bidir.si;
					current_parameters[2] = fits[b].components[k].bidir.l;
					current_parameters[3] = fits[b].components[k].bidir.w;
					current_parameters[4] = fits[b].components[k].bidir.pi;
					current_parameters[5] = fits[b].components[k].bidir.foot_print;
					bootstrapped_parameters[k] 	= current_parameters;
				}
				bootstrapped_parameters 	= sort_bootstrap_parameters(bootstrapped_parameters);
				for (int k = 0 ; k < K; k++){
					double mu_o 	= (segments[i]->parameters[k][0] - segments[i]->start) / scale;
					double si_o 	= segments[i]->parameters[k][1]/scale;
					double l_o 		= scale / segments[i]->parameters[k][2];
					double w_o 		= segments[i]->parameters[k][3];
					double pi_o 	= segments[i]->parameters[k][4];
					double fp_o 	= segments[i]->parameters[k][6]/scale;
					variances[k][0].push_back(abs(mu_o - bootstrapped_parameters[k][0]));
					variances[k][1].push_back(abs(si_o - bootstrapped_parameters[k][1]));
					variances[k][2].push_back(abs(l_o - bootstrapped_parameters[k][2]));
					variances[k][3].push_back(abs(w_o - bootstrapped_parameters[k][3]));
					variances[k][4].push_back(abs(pi_o - bootstrapped_parameters[k][4]));
					variances[k][5].push_back(abs(fp_o - bootstrapped_parameters[k][5]));
					
				}
			}
			for (int k = 0; k < K;k++){
				segments[i]->variances[k] 	= vector<double>(6);
				segments[i]->variances[k][0] 	= get_mean(variances[k][0]);
				segments[i]->variances[k][1] 	= get_mean(variances[k][1]);
				segments[i]->variances[k][2] 	= get_mean(variances[k][2]);
				segments[i]->variances[k][3] 	= get_mean(variances[k][3]);
				segments[i]->variances[k][4] 	= get_mean(variances[k][4]);
				segments[i]->variances[k][5] 	= get_mean(variances[k][5]);
				
			}
		}
	}
}