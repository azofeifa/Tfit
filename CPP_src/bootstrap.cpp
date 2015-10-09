#include "bootstrap.h"
#include "load.h"
#include <vector>
#include "read_in_parameters.h"
#include "model.h"
#include <random>
#include <omp.h>
#include "template_matching.h"
using namespace std;

int sample(double ** CDF, int XN, double sum_N, segment * NS, double pi ){
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
		if (U > pi ){
			s 	= 1;
		}else{
			s 	= 2;
		}
		j 		= 0;
		while (j+1 < XN and CDF[s][j] < U){
			j++;
		}
		NS->X[s][j]++;

		ct++;
	}
	return i;
	
}

void subsample(segment * S, segment * NS ){
	double ** CDF 	= new double*[3];
	NS->X 			= new double*[3];
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
	pi 		= forward_sum / sum_N;
	forward_sum = 0, reverse_sum = 0;
	for (int i = 0; i < S->XN; i++){
		forward_sum+=S->X[1][i];
		reverse_sum+=S->X[2][i];
		CDF[1][i] 	= forward_sum / sum_N;
		CDF[2][i] 	= reverse_sum / sum_N;
	}
	sample( CDF, S->XN, S->N, NS , pi);
	
}

vector<segment *> make_bootstraps(segment * S, params * P){
	int rounds 	= stoi(P->p6["-rounds"]);
	int brounds = stoi(P->p6["-brounds"]);
	int np 		= omp_get_max_threads();
	//make segmnets
	vector<segment *> segments(brounds);

	#pragma omp parallel for num_threads(np)
	for (int r = 0; r < brounds; r++){
		//subsample get new segment

		segment * NS 	= new segment(S->chrom, S->start, S->stop);
		NS->parameters 	= S->parameters;
		subsample(  S,   NS );
		segments[r] 	= NS;
		delete NS;
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


void run_bootstrap_across(vector<segment *> segments, params * P){
	int rounds 		= stoi(P->p6["-rounds"]);
	double scale 	= stod(P->p6["-ns"]);
	int np 			= omp_get_max_threads();
	for (int i = 0; i < segments.size(); i++){
		printf("%d\n",i );
		int K 				= segments[i]->parameters.size();
		vector<double> mu_seeds( K ); 
		double foot_print 	= 0;
		for (int c = 0; c < segments[i]->parameters.size(); c++){
			foot_print+=segments[i]->parameters[c][5];
			mu_seeds[c]=(segments[i]->parameters[c][0] - segments[i]->start) / scale;
		}
		foot_print/=double( K );
		if (segments[i]->N > 0){
			vector<segment *> bootstraps = make_bootstraps(segments[i], P);
			int B 				= bootstraps.size();
			vector<classifier> clfs( B );
			vector<classifier> fits(B);
			#pragma omp parallel for num_threads(np)
			for (int b = 0; b < bootstraps.size(); b++ ){
				double best_ll 	= nINF;
				classifier best_clf;
				segment * S 	= bootstraps[b];
				for (int r = 0 ; r < rounds; r++ ){
					classifier current_clf(K, stod(P->p6["-ct"]), stoi(P->p6["-mi"]), stod(P->p6["-max_noise"]), 
						stod(P->p6["-r_mu"]), stod(P->p6["-ALPHA_0"]), stod(P->p6["-BETA_0"]), stod(P->p6["-ALPHA_1"]), 
						stod(P->p6["-BETA_1"]), stod(P->p6["-ALPHA_2"]) , stod(P->p6["-ALPHA_3"]), false,foot_print );
					current_clf.fit(bootstraps[b], mu_seeds);
					if (r == 0 or current_clf.ll > best_ll){
						best_clf 	= current_clf;
					}
				}
				fits[b]=best_clf;
			}

			for (int b =0 ; b < B; b++){
				vector<vector<double>> bootstrapped_parameters(K);
				for (int k = 0; k < K; k++){
					vector<double> current_parameters(5);
					current_parameters[0] = fits[b].components[k].bidir.mu;
					current_parameters[1] = fits[b].components[k].bidir.l;
					current_parameters[2] = fits[b].components[k].bidir.si;
					current_parameters[3] = fits[b].components[k].bidir.w;
					current_parameters[4] = fits[b].components[k].bidir.pi;
					bootstrapped_parameters[k] 	= current_parameters;
				}
			}


		}
	}
}