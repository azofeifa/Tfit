
#include "load.h"
#include "model.h"
#include <math.h> 
#include <limits>
#include <iostream>
#include <algorithm>
#include "template_matching.h"
#include <fstream>
#include "model_single.h"
#include <random>
#include "omp.h"
using namespace std;

double nINF	=-exp(1000);
double INF 	= exp(1000);


//=========================================================
//Peak Finding

vector<vector<double> > find_peaks(vector<double> values, double ** data){
	vector<vector<double> > peaks;
	for (int i = 1; i < values.size()-1; i++ ){
		if (values[i-1] < values[i] && values[i] > values[i+1]){
			vector<double> current(2); 
			current[0] 	= data[0][i];
			current[1] 	= values[i];
			peaks.push_back(current);
		}
	}
	return peaks; 
}

//=========================================================
//BiDirectional Mappings 

double template_ll(segment * data, int center, 
	int start, int stop, double si, double l ){
	double EMG_ll=0;
	double UNI_ll=0;
	double vl 	= 0.5/(data->X[0][stop] - data->X[0][start]);
	EMG EMG_clf(data->X[0][center], si, l, 1.0, 0.5 );
	for (int i = start; i < stop; i++ ){
		UNI_ll+=(LOG(vl)*data->X[1][i] + LOG(vl)*data->X[2][i]);
		EMG_ll+=(LOG(EMG_clf.pdf(data->X[0][i], 1))*data->X[1][i] + LOG(EMG_clf.pdf(data->X[0][i], -1))*data->X[2][i] );
	}
	return UNI_ll/EMG_ll;
	
}
double window_cov(segment * data, int start, int stop){
	double forward 	= 0;
	double reverse 	= 0; 
	for (int i = start; i < stop; i++ ){
		forward+=data->X[1][i];
		reverse+=data->X[2][i];
	}	
	return LOG(forward)+LOG(reverse) ;
}

vector<vector<double>> bubble_sort(vector<vector<double>> X){ //sort vector of vectors by second
	bool changed=true;
	while (changed){
		changed=false;
		for (int i = 0; i < X.size()-1; i++  )	{
			if (X[i][1] < X[i+1][1]){
				vector<double> copy 	= X[i];
				X[i] 					= X[i+1];
				X[i+1] 					= copy;
				changed=true;
			}
		}
	}
	return X;
}

vector<vector<double>> bubble_sort2(vector<vector<double>> X){ //sort vector of vectors by second
	bool changed=true;
	while (changed){
		changed=false;
		for (int i = 0; i < X.size()-1; i++  )	{
			if (X[i][0] < X[i+1][0]){
				vector<double> copy 	= X[i];
				X[i] 					= X[i+1];
				X[i+1] 					= copy;
				changed=true;
			}
		}
	}
	return X;
}


vector<double> window_search(segment * data, double si, double l, bool kind){
	vector<double> values(data->XN);
	int j,k;
	double window = si + (1.0 / l);
	double vl;
	for (int i = 0; i < data->XN; i++){
		j=i;
		while (j < data->XN && (data->X[0][j] - data->X[0][i]) < window){
			j++;
		}
		k=i;
		while (k >0 and (data->X[0][i] - data->X[0][k]) < window){
			k--;
		}
		if (kind){
			vl =	template_ll(data, i, k,j, si, l);
		}else{
			vl =	window_cov(data, k, j);			
		}
		values[i] 	= vl;
	}

	return values;
}
vector<double> take_avg(vector<double> X1, vector<double> X2 ){
	vector<double> avg;
	if (X1.size() != X2.size()){
		cout<<"problem with template matchers...."<<endl;
		return avg;
	}
	for (int i = 0; i < X1.size(); i++){
		avg.push_back((X1[i] + X2[i])/2.);
	}
	return avg;
}


vector<double> peak_bidirs(segment * data){
	vector<double> coverage_values 	= (window_search(data, 1, 0.1, 0));
	vector<double> template_values 	= (window_search(data, 1, 0.1, 1));
	vector<double> average_values 	= take_avg(coverage_values, template_values);
	
	vector<vector<double>> peaks 	= find_peaks(average_values, data->X);
	vector<double> centers ;
	if (peaks.size()==0){
		return centers;
	}
	peaks 							= bubble_sort(peaks);
	for (int i = 0 ; i < peaks.size(); i++)	{
		centers.push_back(peaks[i][0]);
	}
	return centers;
}

int sample_centers(vector<double> centers, double p){
	random_device rd;
	mt19937 mt(rd());
	default_random_engine generator;
	uniform_int_distribution<int> distribution(0,centers.size()-1);
	int i 	= distribution(mt);
	//int s 	= centers.size()-1;
	// if (i > s){
	// 	return s;
	// }
	return i;
}

void sample_variance(double ** X, double mu, double * array, int j,int k, double foot_print){
	double N_forward, N_reverse;
	array[0] 	= 0;
	array[1] 	= 0;
	N_forward=0, N_reverse=0;
	
	for (int i =k; i < j;i++ ){
		array[0]+=pow((X[0][i]-foot_print)-mu,2)*X[1][i];
		array[1]+=pow((X[0][i]+foot_print)-mu,2)*X[2][i];
		N_forward+=X[1][i];
		N_reverse+=X[2][i];
	}
	double s;
	array[0] 	= array[0] / N_forward;
	array[1] 	= array[1] / N_reverse;
}
double sample_variance_single(double ** X, double mu, int j, int k){
	double N 	= 0;
	double S 	= 0;
	for (int i =k; i < j; i++){
		S+=pow(mu - X[0][i], 2)*X[1][i];
		N+=X[1][i];
	}
	return S/N;
}

void sample_variance2(double ** X, double mu, double * array, int j,int k){
	double N_forward, N_reverse;
	array[0] 	= 0;
	array[1] 	= 0;
	N_forward=0, N_reverse=0;
	for (int i =k; i < j;i++ ){
		if (X[0][i] > mu){
			array[0]+=pow(X[0][i]-mu,2)*X[1][i];
			N_forward+=X[1][i];
		}
		if (X[0][i]<mu){
			array[1]+=pow(X[0][i]-mu,2)*X[2][i];
			N_reverse+=X[2][i];
	
		}
	}
	double s;
	array[0] 	= array[0] / N_forward;
	array[1] 	= array[1] / N_reverse;	
}
void sample_skew(double **X , double mu, double * sample_variances, double * sample_skews, int j, int k, double foot_print){
	double N_forward, N_reverse;
	sample_skews[0] 	= 0;
	sample_skews[1] 	= 0;
	for (int i = k; i<j;i++){
		sample_skews[0]+=pow((X[0][i]-foot_print)-mu,3)*X[1][i];
		sample_skews[1]+=pow((X[0][i]+foot_print)-mu,3)*X[2][i];
		N_forward+=X[1][i];
		N_reverse+=X[2][i];
	}
	sample_skews[0]/=N_forward;
	sample_skews[1]/=N_reverse;
	sample_skews[0]/=pow(sample_variances[0], 3/2.);
	sample_skews[1]/=pow(sample_variances[1], 3/2.);	
}
double moment_estimate_lambda(double ** X, double mu,int j, int k, double foot_print ){
	double N_forward, N_reverse;
	double means[2];
	means[0]=0;
	means[1]=0;
	for (int i = k; i < j ;i++){
		means[0]+=(X[0][i]-foot_print)*X[1][i];
		N_forward+=X[1][i];
		means[1]+=(X[0][i]+foot_print)*X[2][i];
		N_reverse+=X[2][i];		
	}
	return 1.0 / (0.5*((means[0]/N_forward) - (means[1]/N_reverse) ));
}
double moment_estimate_sigma(double * sample_variances, double l){
	double std 	= 0.5*(sqrt(sample_variances[0]) + sqrt(sample_variances[1]));
	return std 	= std - (1. / l);
}

double get_mean(double ** X, int j,int k,int s, double foot_print){
	double SUM,N;
	SUM=0;N=0;
	int st;
	if (s == 1){
		st 	= s;
	}else{
		st 	= -1;
	}
	for (int i = k ; i<j; i++){
		SUM+=(X[0][i]-st*foot_print)*X[s][i];
		N+=X[s][i];
	}
	return SUM/N;
}

double BIC(double ** X,  double * avgLL, double * variances,double * lambdas, 
	double ** skews, double mu, int i, int j, int k , int single, double fp_res, double * densities_f, double * densities_r, double scale){
	double score 	= 0;
	double a 		= X[0][k];
	double b 		= X[0][j];	
	double fp_a 	= 0;
	double fp_b 	= 1000;
	double fp_delta = (fp_b - fp_a) / fp_res;
	double foot_print;
	if (not single){
		double emg_ll, uni_ll;
		emg_ll=0, uni_ll=0;
		variances[i] 	= 0;
		lambdas[i] 		= 0;
		avgLL[i] = nINF;
		
		double * best_skew 					= new double[2];
		best_skew[0]=0, best_skew[1]=0;
		skews[i][0]  	= best_skew[0], skews[i][1]= best_skew[1];
	
		double arg_l, arg_si, arg_bic, arg_ll; 
		arg_l=0, arg_si=0;
		arg_ll=nINF, arg_bic=INF;
		double N_pos=0;
		double N_neg=0;
		double N=0;
		for (int i = k; i < j; i++ ){
			N_pos+=X[1][i];
			N_neg+=X[2][i];
		}
		densities_f[i] 	= N_pos / (b-a);
		densities_r[i] 	= N_neg / (b-a);
		if (N_pos == 0 or N_neg ==0){
			return 0;
		}
		N 	  = (N_pos+N_neg);
		double pi 	= N_pos / N;
		uni_ll= LOG(pi/ (b-a) )*N_pos + LOG((1-pi)/(b-a))*N_neg;
		for (int fp =0; fp< fp_res; fp++){
			foot_print 	= (fp_a + fp_delta*fp)/scale;
			//(get_mean(X, j,k,1, foot_print) + get_mean(X, j,k,2, foot_print))*0.5  ;
			double * current_sample_variances 	= new double[2];
			double * current_sample_skews 		= new double[2];
			sample_variance(X, mu, current_sample_variances, j,k, foot_print);
			sample_skew(X, mu, current_sample_variances, current_sample_skews, j,k,foot_print);
			double l 	= moment_estimate_lambda(X,mu,j,k, foot_print);
			double sii 	= moment_estimate_sigma(current_sample_variances, l);
			if (l > 0 and sii > 0){
				EMG EMG_clf(mu, sii, l, 1.0, 0.5 );
				emg_ll=0;
				for (int i = k; i < j; i++ ){
					emg_ll+=(LOG(EMG_clf.pdf((X[0][i]- foot_print),1))*X[1][i] + LOG(EMG_clf.pdf((X[0][i]+foot_print),-1))*X[2][i]);	
				}
				if (emg_ll > arg_ll){
					arg_ll=emg_ll, arg_bic= (-2*uni_ll + 1*LOG(N) ) / (-2*emg_ll + 3*LOG(N));
					arg_l=l, arg_si=sii;
					best_skew 	= current_sample_skews;
				}	
			}
		}
		if (arg_l==0 or arg_si==0 or arg_ll==nINF){

			return 0;
		}

		variances[i] 	= arg_si;
		lambdas[i] 		= arg_l;
		avgLL[i] 		= arg_ll / N;
		skews[i][0]  	= best_skew[0], skews[i][1]= best_skew[1];
		score 			= arg_bic;
	}else{
		lambdas[i] 		= 0;
		variances[i] 	= 0;
		skews[i][0] = 0,skews[i][1] = 0;
		double MU 		= mu;
		double SI 		= sample_variance_single(X,MU, j,k);
		variances[i] 	= SI;
		
		NORM clf(MU, sqrt(SI), 1);
		double n_ll, uni_ll, N;
		n_ll=0, uni_ll=0, N=0;
		for (int i = k; i < j; i++ ){
			n_ll+=(LOG(clf.pdf(X[0][i] ))*X[1][i] );
			uni_ll+=(LOG(	1.0 / (b-a) )*X[1][i]  ); 
			N+=(X[1][i] );
		}
		double BIC_score_emg 	= -2*n_ll + 5*LOG(N) ;
		double BIC_score_uni 	= -2*uni_ll + 1*LOG(N) ;
		avgLL[i] 				= n_ll / N;
		score 	=  BIC_score_uni / BIC_score_emg;
	}
	return score;


	

}

void BIC_template(segment * data, double * avgLL, double * BIC_values, double * densities, double * densities_r,
	double * variances,double * lambdas, double ** skews ,double window, int np, int single,double foot_res,double scale){
	double vl;
	int NN 	= int(data->XN);
	int threads  	= omp_get_max_threads();
	#pragma omp parallel for num_threads(threads)
	for (int i = 0; i < NN; i++){
		int j=i;
		while (j < data->XN && (data->X[0][j] - data->X[0][i]) < window){
			j++;
		}
		int k=i;
		while (k >0 and (data->X[0][i] - data->X[0][k]) < window){
			k--;
		}
		if (k > 0 and j < data->XN ){
			if (not single){
				BIC_values[i] 	= BIC(data->X, avgLL, variances, lambdas, skews, data->X[0][i], i, j, k,  single, foot_res, densities, densities_r,scale );
			}else{
				BIC_values[i] 	= BIC(data->X, avgLL, variances, lambdas, skews, data->X[0][i], i, j, k,  single, foot_res, densities, densities_r,scale );	
			}
				
		}else{
			BIC_values[i] 	= 0;
		}
	}
}




void run_global_template_matching(vector<segment*> segments, 
	string out_dir,  double res, double density,
	double scale, double ct, int np, double skew, int single, ofstream& log_file){
	
	ofstream FHW;
	ofstream FHW_intervals;
	
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
	if (not single){
		window_a 	= 500;
		window_b 	= 3000;
	}else{
		window_a 	= 2000;
		window_b 	= 4000;
	
	}
	int all 	= 0;
	double fp_res = 10;
	double window_delta = (window_b-window_a)/res;
	window 		= 1500/scale;
	log_file<<"(template_matching) beginning to run template matching...";
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
			double current_BIC_score 	= INF;
			vector<double> arg_bic;
			for (int i = 0; i < intervals.size(); i++){
				if (intervals[i][2] < current_BIC_score){
					arg_bic 	= intervals[i];
				}
			}
			if (arg_bic.empty()){
				printf("WHAT?\n");
			}
			return arg_bic;
		}
		
	};
	int mj;
	vector<merged> mergees;
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
		for (int w = 0 ; w<res; w++){
			window 		= (window_a+window_delta*w)/scale;
			BIC_template(segments[i], avgLL, BIC_values, densities, densities_r, variances, lambdas,skews, window, np, single, fp_res,scale);
			mj 	= 0;
			//write out contigous regions of up?
			for (int j = 1; j<segments[i]->XN-1; j++){
				if (avgLL[j-1]< avgLL[j] and avgLL[j] > avgLL[j+1]){
					if (BIC_values[j] >=ct  and densities[j]>window*0.75  and densities_r[j]>window*0.75 and skews[j][0] >= skew and skews[j][1] <= -skew){
						if (lambdas[j]>0){
							start 		= int(segments[i]->X[0][j]*scale+segments[i]->start - ((variances[j]/2.)+(1.0/lambdas[j]))*scale);
							stop 		= int(segments[i]->X[0][j]*scale+segments[i]->start + ((variances[j]/2.)+(1.0/lambdas[j]))*scale);
						}else{
							start 		= int(segments[i]->X[0][j]*scale+segments[i]->start - (variances[j]/2.)*scale);
							stop 		= int(segments[i]->X[0][j]*scale+segments[i]->start + (variances[j]/2.)*scale);	
						}
						current[0] 	= double(start), current[1]=double(stop), current[2]=avgLL[j], current[3]=(variances[j]/4.)*scale, current[4]=(2/lambdas[j])*scale;
						
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
		
			
		}
		
		scores.clear();
		for (int m = 0; m < mergees.size(); m++){
			scores.push_back(mergees[m].get_best());
		}
		all+=int(scores.size());
		for (int j = 0; j < scores.size();j++){

			center 		= int((scores[j][1]+scores[j][0]) / 2.);
			//want to insert into
			vector<double> bounds(2);
			bounds[0] 	= scores[j][0], bounds[1]=scores[j][1];
			segments[i]->bidirectional_bounds.push_back(bounds);	
		}
		scores.clear();
	}
	log_file<<"done, found: " + to_string(all) + " potential sites of bidirectional transcription\n";
}


void optimize(map<string, interval_tree *> I, 
	vector<segment*> segments, double scale, int res, string out_dir, string spec_chrom, int np){
	double window_a 		= 500/scale;
	double window_b 		= 3000/scale;
	double window_delta 	= (window_b - window_a) / res;
	
	double density_a 		= (100/scale);
	double density_b 		= (5000/scale);
	double density_delta 	= (density_b - density_a) / res;
	double skew_a 			= -1;
	double skew_b 			= 1.;
	double skew_delta 		= (skew_b-skew_a)/res;
	
	double ct_a 			= 0.9;
	double ct_b 			= 1.2;
	double ct_delta 		= (ct_b - ct_a) / res;
	
	double accuracy[res+1][res+1][res+1];
	
	double density, ct, peak, skew, window;

	int prev, prev_start, current, stop;
	int bidir_hits 	= 0;
	int all_bidirs 	= 0;
	int no_hits 	= 0;
	int min_start 	= 0;
	int max_stop 	= 0;
	double F1 		= 0;
	int arg_j, arg_k, arg_l, arg_v;
	bool HIT;
	bool FOUND=false;
	double skew_sum_f, skew_N_f;
	double skew_sum_r, skew_N_r;

	bool GOOD=true;
	printf("---------------------\n");
	printf("performing optimization\n");
	for (int i = 0; i < segments.size(); i++){ //these segments are indexed by chromosome
		if (segments[i]->chrom==spec_chrom){
			FOUND 	= true;
			for (int v =0 ; v < res; v++){
				window 	= window_a + window_delta*v;
				double * avgLL 		= new double[int(segments[i]->XN)];
		
				double * BIC_values = new double[int(segments[i]->XN)];
				double * densities= new double[int(segments[i]->XN)];
				double * densities_r= new double[int(segments[i]->XN)];
				
				double * variances= new double[int(segments[i]->XN)];
				double * lambdas= new double[int(segments[i]->XN)];
				double ** skews = new double*[int(segments[i]->XN)] ;
				for (int t =0; t < segments[i]->XN; t++ ){
					skews[t] 	= new double[2];
				}
				BIC_template(segments[i], avgLL, BIC_values, densities, densities_r, variances, lambdas, skews, window,np,0,0,scale);
				
				for (int k =0; k <= res; k++){
					density 	= density_a+density_delta*(k);
					for (int l =0; l <= res; l++){
						bidir_hits 	= 0;
						all_bidirs 	= 0;
						min_start 	= segments[i]->X[0][0]*scale +segments[i]->start ;
			
						prev=-1, current=0;

						ct 	= ct_a 	+ ct_delta*(l);
						for (int j = 1; j<segments[i]->XN-1; j++){
							if (avgLL[j-1]< avgLL[j] and avgLL[j] > avgLL[j+1]){
								if (BIC_values[j] >= ct and densities[j]>density and densities_r[j] > density 
									and skews[j][0] > 0.1 and skews[j][1]< -0.1  ){
									prev_start 	= segments[i]->X[0][j]*scale+segments[i]->start - ((variances[j]/2.)+(1./lambdas[j]))*scale;
									stop 		= segments[i]->X[0][j]*scale+segments[i]->start + ((variances[j]/2.)+(1./lambdas[j]))*scale;
									HIT 		= I[segments[i]->chrom]->find(prev_start, stop);
									if (HIT){
										bidir_hits++;
									}
									all_bidirs++;
								}
							}
							max_stop 	= segments[i]->X[0][j]*scale +segments[i]->start;
						}
						int single_hits 	= I[segments[i]->chrom]->get_hits(true, min_start, max_stop);
						int all_hits 		= I[segments[i]->chrom]->get_hits(false, min_start, max_stop);
						int peaks_total 	= I[segments[i]->chrom]->get_total(min_start, max_stop);
						double precision 	= double(bidir_hits)/double(all_bidirs);
						double recall 		= double(bidir_hits)/double(peaks_total);
						if ( 2*((precision*recall)/(precision+recall) ) > F1  ){
							F1			= 2*((precision*recall)/(precision+recall) ) ;
							arg_l=l, arg_k=k , arg_v=v;
							GOOD 		= true;
						}
						I[segments[0]->chrom]->reset_hits();
					}	
				printf("Best F1 (so far): %.3g, window: %.3gKB, density: %.3gKB, ct: %.3g, percent done: %d  \r", F1, 
				( (window_a + window_delta*arg_v)*scale/1000.  ), (scale*(density_a+density_delta*(arg_k))/1000.), (ct_a 	+ ct_delta*(arg_l)),
				int(100*double(v+1)/(res+1)) );
				cout<<flush;
				}
			}
			
		}
	}
	if ( not GOOD ){
		printf("optimization failed\n");
	}else{
		printf("\nFinished Optimization");
		printf("\n---------------------\n");
		

	}
		
}










