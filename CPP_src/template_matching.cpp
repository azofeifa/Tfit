
#include "load.h"
#include "model.h"
#ifdef USING_ICC
#include <mathimf.h>
#else
#include <math.h>   
#endif
#include <limits>
#include <iostream>
#include <algorithm>
#include "template_matching.h"
#include <fstream>
#include <random>
#include "omp.h"
#ifdef USING_ICC
#include <aligned_new>
#endif
using namespace std;

double nINF	=-exp(1000);
double INF 	= exp(1000);


vector<vector<double>> bubble_sort3(vector<vector<double>> X){ //sort vector of vectors by second
	bool changed=true;
	if (X.size()<2){
		return X;
	}
	while (changed){
		changed=false;
		for (int i = 0; i < X.size()-1; i++  )	{
			if (X[i][2] < X[i+1][2]){
				vector<double> copy 	= X[i];
				X[i] 					= X[i+1];
				X[i+1] 					= copy;
				changed=true;
			}
		}
	}
	return X;
}

int sample_centers(vector<double> centers, double p){
	random_device rd;
	mt19937 mt(rd());
	default_random_engine generator;
	uniform_int_distribution<int> distribution(0,centers.size()-1);
	int i 	= distribution(mt);
	return i;//required in model.o (ugh...)
}

double BIC2(double ** X,  double * avgLL, double * variances,double * lambdas, 
	double ** skews, double N_pos, double N_neg, double S_pos, 
		double S_neg, double S2_pos, double S2_neg,double mu, int j,int k,int i, double scale ){
	double argBIC 	= 0;
	double arg_si 	= 0;
	double arg_l 	= 0;
	double arg_ll 	= nINF;
	double fp_a 	= 0;
	double fp_b 	= 500/scale;
	double fp_res 	= 5;
	double fp_delta = (fp_b-fp_a) /fp_res;
	double N 	= N_neg + N_pos;
	double pi 	= (N_pos) / (N_neg + N_pos);
	double pi2 	= (N_pos+10000) / (N_neg + N_pos+20000);
	double a 	= X[0][j], b=X[0][k];
	double uni_ll= LOG(pi/ pow(b-a,1) )*N_pos + LOG((1-pi)/pow(b-a,1))*N_neg;
	for (int fp = 0; fp < fp_res; fp++){

		double foot_print 	= fp*fp_delta;
		
		double l  	= 1./ (0.5*((S_pos / N_pos) - (S_neg / N_neg)) - foot_print);
		double sv_f = sqrt((S2_pos - (2*(mu )*S_pos) + (N_pos*pow((mu ),2)))/N_pos);
		double sv_r = sqrt((S2_neg - (2*(mu )*S_neg) + (N_neg*pow((mu ),2)))/N_neg);
		double si 	= 0.5*(sv_f + sv_r) - (1. / l);
		if (l > 0 and si > 0){

			EMG EMG_clf(mu, si, l, 1.0, pi );
			double emg_ll=0;
			for (int i = j; i < k; i++ ){
				emg_ll+=(LOG(EMG_clf.pdf((X[0][i] ),1))*X[1][i] + LOG(EMG_clf.pdf((X[0][i] ),-1))*X[2][i]);	
			}	
			double currBIC= (-2*uni_ll + 1*LOG(N) ) / (-2*emg_ll + 3*LOG(N));
			if (currBIC > argBIC){
				argBIC=currBIC, arg_si=si, arg_l=l, arg_ll=emg_ll;
			}
		}
	}
	variances[i] 	= arg_si;
	lambdas[i] 		= arg_l;
	avgLL[i] 		= arg_ll / N;
	skews[i][0]  	= 0, skews[i][1]= 0;

	return argBIC;
}

double get_ll(double ** X, double mu, double w, double pi, double l, int j, int k){
	double emg_ll 	= 0;
	EMG EMG_clf(mu, 1.0, 0.1, w, pi  );
	for (int i = j; i < k;i++ ){
		emg_ll+=LOG(EMG_clf.pdf(X[0][i],1) + (1.0-w)*pi*(1.0/l) )*X[1][i] + LOG(EMG_clf.pdf(X[0][i],-1) + (1.0-w)*(1.0-pi)*(1.0/l) )*X[2][i];
	}
	return emg_ll;

}


double BIC3(double ** X, int j, int k, int i,
	double N_pos, double N_neg, double * avgLL, double * variances,double * lambdas, 
	double ** skews){

	double w = 0.01;
	double delta 	= 0.25;
	double N 		= N_pos + N_neg;
	double pi 	= (N_pos+10000) / (N_neg + N_pos+20000);
	double a 	= X[0][j], b=X[0][k];
	double uni_ll= LOG(pi/  (b-a ) )*N_pos + LOG((1-pi)/ (b-a ))*N_neg;
	double best_emg_ll 	= nINF;
	double 		l = b-a;
	double best_w 	= 0.0;
	w 				= 0.95;
	double fp 		= 0.0;

	double emg_ll 	= 0;
	EMG EMG_clf(X[0][i], 1.0, 0.1, w, pi  );
	EMG_clf.foot_print 	= 1.0;

	for (int i = j; i < k;i++ ){
		emg_ll+=LOG(EMG_clf.pdf(X[0][i],1) + (1.0-w)*pi*(1.0/l) )*X[1][i] + LOG(EMG_clf.pdf(X[0][i],-1) + (1.0-w)*(1.0-pi)*(1.0/l) )*X[2][i];
	}

	
	double emg_ratio 			= (-2*uni_ll + LOG(N)) / (-2*emg_ll + 1*LOG(N));

	variances[i] 	= 1.0;
	lambdas[i] 		= 1.0;
	avgLL[i] 		= best_emg_ll / N;
	skews[i][0]  	= 0, skews[i][1]= 0;
	return emg_ratio;



}













void BIC_template(segment * data, double * avgLL, double * BIC_values, double * densities, double * densities_r,
	double * variances,double * lambdas, double ** skews ,double window, int np, int single,double foot_res,double scale){
	double vl;
	int NN 	= int(data->XN);
	int threads  	= omp_get_max_threads();
	int counts 		= NN / threads;
	#pragma omp parallel num_threads(threads)
	{
		int tid 	= omp_get_thread_num();
		int start 	= tid*counts;
		int stop 	= (tid+1)*(counts);
		if (tid+1 == threads){
			stop 	= NN;
		}
		int j = start, k =start;
		double N_pos=0,N_neg=0;
		double S_pos=0, S_neg=0;
		double S2_pos=0, S2_neg=0;
		double total_density;
		for (int i = start; i < stop; i++){
			while (j < data->XN and (data->X[0][j] - data->X[0][i]) < -window){
				N_pos-=data->X[1][j];
				N_neg-=data->X[2][j];
				S_pos-=(data->X[0][j]*data->X[1][j]);
				S_neg-=(data->X[0][j]*data->X[2][j]);
				
				S2_pos-=(pow(data->X[0][j],2)*data->X[1][j]);
				S2_neg-=(pow(data->X[0][j],2)*data->X[2][j]);
				
				j++;
			}
			while (k < data->XN and (data->X[0][k] - data->X[0][i]) < window){
				N_pos+=data->X[1][k];
				N_neg+=data->X[2][k];

				S_pos+=(data->X[0][k]*data->X[1][k]);
				S_neg+=(data->X[0][k]*data->X[2][k]);

				S2_pos+=(pow(data->X[0][k],2)*data->X[1][k]);
				S2_neg+=(pow(data->X[0][k],2)*data->X[2][k]);
				k++;
			}
			if (k < data->XN  and j < data->XN and k!=j and N_neg > 0 and N_pos > 0 ){
				total_density 	= (N_pos / (data->X[0][k] - data->X[0][j])) + (N_neg / (data->X[0][k] - data->X[0][j]));
				densities[i] 	= N_pos ;
				densities_r[i] 	= N_neg ;
				
				// BIC_values[i] 	= BIC(data->X, avgLL, variances, 
				// 	lambdas, skews, data->X[0][i], i, k, j,  single, foot_res, 
				// 	densities, densities_r,scale , window, N_pos, N_neg);	
				double mu 		= (data->X[0][k] + data->X[0][j]) /2.;
				//double CENTER 	= (mu*100 + data->start);


				// BIC_values[i] 	=  BIC2(data->X, avgLL, variances, 
				// 	lambdas, skews, N_pos,  N_neg, S_pos, 
			 // 		S_neg, S2_pos, S2_neg, mu, j, k,i, scale );
				BIC_values[i] 	= BIC3(data->X,  j,  k,  i, N_pos,  N_neg,avgLL, variances, lambdas, skews);




				double II 		= data->X[0][i]*scale + data->start;
				// if   (CENTER < 1144764 and CENTER > 1139533 and data->chrom == "chr1"){
				// 	printf("%f,%f,%f\n", BIC_values[i],CENTER, densities[i],densities_r[i] );
				// }
				
			
					
			}else{
				BIC_values[i] 	= 0;
				densities[i] 	= 0;
				densities_r[i] 	= 0;
			}
		}
	}
}

void run_global_template_matching(vector<segment*> segments, 
	string out_dir,  double res, double density,
	double scale, double ct, int np, double skew, int single){
	
	ofstream FHW;
	ofstream FHW_intervals;
	ofstream FHW_scores;
	string annotation;
	int prev, prev_start, stop;
	int N;
	int start, center;
	vector<vector<double> > scores;
	vector<double> current(5);
	vector<string> INFOS;
	double window, foot_print;
	double window_a;
	double window_b;
	window_a 	= 1000;
	window_b 	= 1500;

	int all 	= 0;
	double fp_res = 10;
	double window_delta = (window_b-window_a)/res;
	window 		= 2500/scale;
	//now we want to merge all of these overlaps...
	struct merged{
	public:
		double start, stop;
		vector<vector<double>> intervals;
		void add(vector<double> insertee){
			start 		= min(start, insertee[0]), stop= max(stop, insertee[1]);
			intervals.push_back(insertee);
		}
		merged(){};
		merged(vector<double> insertee ){
			start=insertee[0], stop=insertee[1];
			intervals.push_back(insertee);
		};
		vector<double> get_best(){
			double current_BIC_score 	= 0;
			vector<double> arg_bic 		= intervals[0];
			for (int i = 0; i < intervals.size(); i++){
				if (intervals[i][2] > current_BIC_score){
					arg_bic 	= intervals[i];
				}
			}
			if (arg_bic.empty()){
				printf("%d\n",intervals.size() );
				printf("WHAT?\n");
			}
			return arg_bic;
		}
		
	};
	int mj;
	int threads  	= omp_get_max_threads();
	vector<merged> mergees;
	//#pragma omp parallel for num_threads(threads)
	for (int i = 0; i < segments.size(); i++){
		double * avgLL 			= new double[int(segments[i]->XN)];
		double * BIC_values 	= new double[int(segments[i]->XN)];
		double * densities 		= new double[int(segments[i]->XN)];
		double * densities_r 	= new double[int(segments[i]->XN)];
		double * variances 		= new double[int(segments[i]->XN)];
		double * lambdas 		= new double[int(segments[i]->XN)];
		double ** skews 		= new double*[int(segments[i]->XN)] ;
		for (int t =0; t < segments[i]->XN; t++ ){
			skews[t] 	= new double[2];
		}
		mergees.clear();


		for (int w = 0 ; w<1; w++){
			double l 		=  segments[i]->XN*scale;
			window 			= (window_a+window_delta*w)/(scale);
			window 			= 1000/scale;
			double ef 		= segments[i]->fN*( 2*(window_a+window_delta*w)  /l);
			double er 		= segments[i]->rN*( 2*(window_a+window_delta*w) /l);
			double stdf 	= sqrt(ef*(1- (  2*(window_a+window_delta*w) /l  ) )  );
			double stdr 	= sqrt(er*(1- (  2*(window_a+window_delta*w) /l  ) )  );
//			printf("%f,%f,%f,%f\n",ef,er, stdf,stdr );
			BIC_template(segments[i], avgLL, BIC_values, densities, densities_r, variances, lambdas,skews, window, np, single, fp_res,scale);

			mj 	= 0;
			//write out contigous regions of up?
			for (int j = 1; j<segments[i]->XN-1; j++){
				

				if (BIC_values[j] >=ct and densities[j] > ef + 1*stdf  and densities_r[j]> er + 1*stdr    ){
					// start 		= int(segments[i]->X[0][j]*scale+segments[i]->start - ((variances[j]/2.)+(1.0/lambdas[j]))*scale);
					// stop 		= int(segments[i]->X[0][j]*scale+segments[i]->start + ((variances[j]/2.)+(1.0/lambdas[j]))*scale);
					

					start 		= int(segments[i]->X[0][j]*scale+segments[i]->start - window*scale*0.5);
					stop 		= int(segments[i]->X[0][j]*scale+segments[i]->start + window*scale*0.5);
					
					current[0] 	= double(start), current[1]=double(stop), current[2]=BIC_values[j], current[3]=(variances[j]/4.)*scale, current[4]=(2/lambdas[j])*scale;
					
					merged M(current);
					int N 		= mergees.size();
					mj=0;
					while (mj < N and mergees[mj].stop < M.start){
						mj++;
					}
					if (mj < N and mergees[mj].stop > M.start and mergees[mj].start < M.stop){
						mergees[mj].add(current);
					}else{
						mergees.insert(mergees.begin() + mj, M);
						
					}
				}
			
			}
		
			
		}
		
		scores.clear();
		for (int m = 0; m < mergees.size(); m++){
			vector<double> info(3);
			vector<double> BEST 	= mergees[m].get_best();
			info[0] 	= mergees[m].start, info[1] = mergees[m].stop, info[2] = BEST[2];
			scores.push_back(info);

		}
		all+=int(scores.size());
		for (int j = 0; j < scores.size();j++){

			center 		= int((scores[j][1]+scores[j][0]) / 2.);
			//want to insert into
			vector<double> bounds(3);
			bounds[0] 	= scores[j][0], bounds[1]=scores[j][1], bounds[2]=scores[j][2];
			segments[i]->bidirectional_bounds.push_back(bounds);	
		}
		//sort
		segments[i]->bidirectional_bounds 	= bubble_sort3(segments[i]->bidirectional_bounds);


		scores.clear();
	}
}


void noise_global_template_matching(vector<segment*> segments, double scale){

	double window, foot_print;
	double window_a;
	double window_b;
	double res 	= 2;
	
	double window_delta = (window_b-window_a)/res;
	
	int np  		= omp_get_max_threads();
	int ct 			= 0;
	int ct_thresh 	= 200;
	random_device rd;
	mt19937 mt(rd());
	uniform_real_distribution<double> sample( 0,1 );
	
	for (int i = 0; i < segments.size(); i++){
		double * avgLL 			= new double[int(segments[i]->XN)];
		double * BIC_values 	= new double[int(segments[i]->XN)];
		double * densities 		= new double[int(segments[i]->XN)];
		double * densities_r 	= new double[int(segments[i]->XN)];
		double * variances 		= new double[int(segments[i]->XN)];
		double * lambdas 		= new double[int(segments[i]->XN)];
		double ** skews 		= new double*[int(segments[i]->XN)] ;
		for (int t =0; t < segments[i]->XN; t++ ){
			skews[t] 	= new double[2];
		}
		double l 		=  segments[i]->XN*scale;
		window 			= 1000/(scale);
		
		double ef 		= segments[i]->fN*( 2*(1000)  /l);
		double er 		= segments[i]->rN*( 2*(1000) /l);
		double stdf 	= sqrt(ef*(1- (  2*(1000) /l  ) )  );
		double stdr 	= sqrt(er*(1- (  2*(1000) /l  ) )  );
		BIC_template(segments[i], avgLL, BIC_values, densities, densities_r, 
			variances, lambdas,skews, window, np, 0, 5,scale);
		bool f 	= false;
		int start=0, stop=0;

		for (int j = 1; j<segments[i]->XN-1; j++){
			if ( (BIC_values[j] > 0.05 or densities[j] > ef or densities_r[j] > er)  and  f ){
				stop 	= segments[i]->X[0][j]*scale + segments[i]->start;
				if (stop - start > 1000 and sample(mt) < 0.05 and ct < ct_thresh ){
					vector<double> current(3);
					current[0]=double(start), current[1]=double(stop), current[2]=BIC_values[j-1];
					segments[i]->bidirectional_bounds.push_back(current);		
					ct+=1;
				}
				f 		= false;	
			}
			else if (BIC_values[j] < 0.05 and densities[j] < ef and densities_r[j] < er and !f){
				f 		= true;
				start 	= segments[i]->X[0][j]*scale + segments[i]->start;
			}
		}
		segments[i]->bidirectional_bounds 	= bubble_sort3(segments[i]->bidirectional_bounds);

	}
	
}












