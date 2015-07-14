#include "load.h"
#include "model.h"
#include "template_matching.h"
#include <math.h> 
#include <limits>
#include <iostream>
#include <algorithm>
#include <unistd.h>
#include <random>
//=============================================
//Helper functions
double IN(double x){ //Standard Normal PDF 
	return exp(-pow(x,2)*0.5)/sqrt(2*M_PI);
	
}
double IC(double x){ //Standard Normal CDF
	return 0.5*(1+erf(x/sqrt(2)));
}


double R(double x){ //Mills Ratio
	if (x > 8){
		return 1.0 / x;
	}
	double N = IC(x);
	double D = IN(x);
	if (D < pow(10,-15)){ //machine epsilon
		return 1.0 / pow(10,-15);
	}
	return exp(log(1. - N)-log(D));
}

bool checkNumber(double x){
	if (isfinite(x)){
		return true;
	}
	return false;

}
//=============================================
// Uniform Noise Class
double NOISE::pdf(double x, int strand){
	if (strand == 1){
		return (w*pi) / (b-a);
	}
	return (w*(1-pi)) / (b-a);
}
NOISE::NOISE(){}
NOISE::NOISE(double A, double B, double W, double PI){
	a=A;
	b=B;
	w=W;
	pi=PI;
}

//=============================================
//Uniform Class
double UNI::pdf(double x, int strand){ //probability density function
	double p;
	if (strand != st){
		return 0.;
	}
	if ((a+delta_a)<= x and x<=(b + delta_b)){
		p=w / abs((b+delta_b)-(a+delta_a));
		if (checkNumber(p)){
			return p;
		}
		return 0;
	}
	return 0;
}

string UNI::print(){
	string text = ("U: " + to_string(a) + "," + to_string(b) 
	+ "," + to_string(w) + "," + to_string(pi));
	return text;
}	
UNI::UNI(double start, double stop, double w_i, int strand){

	a=start;
	b=stop;
	w=w_i;
	st=strand;
	if (st==1){
		pi=1.;
	}else{
		pi=0;
	}
	delta_a=0;
	delta_b=0;
	ri_forward=0, ri_reverse=0;

}

UNI::UNI(){} //empty constructor


//=============================================
//Exponentially Modified Gaussian class

string EMG::print(){
	string text 	= ("N: " + to_string(mu)+","+to_string(si)
		+ "," + to_string(l) + "," + to_string(w) + "," + to_string(pi));
	return text;
}

EMG::EMG(){}//empty constructor

EMG::EMG(double MU, double SI, double L, double W, double PI ){
	mu 	= MU;
	si 	= SI;
	l  	= L;
	w 	= W;
	pi 	= PI;
}


double EMG::pdf(double z, int s ){
	double vl 		= (l/2)*(s*2*(mu-z) + l*pow(si,2));
	double p;
	if (vl > 100){ //potential for overflow, inaccuracies
		p 			= l*IN((z-mu)/si)*R(l*si - s*((z-mu)/si));
	}else{
		p 			= (l/2)*exp(vl)*erfc((s*(mu-z) + l*pow(si ,2))/(sqrt(2)*si));
	}
	vl 				= p*w*pow(pi, max(0, s) )*pow(1.-pi, max(0, -s) );
	if (checkNumber(vl)){
		return vl;
	}
	return 0.;

}
double EMG::EY(double z, int s){
	return max(0. , s*(z-mu) - l*pow(si, 2) + (si / R(l*si - s*((z-mu)/si))));
}
double EMG::EY2(double z, int s){

	return pow(l,2)*pow(si,4) + pow(si, 2)*(2*l*s*(mu-z)+1 ) + pow(mu-z,2) - ((si*(l*pow(si,2) + s*(mu-z)))/R(l*si - s*((z-mu)/si) ));
	
}

//=========================================================
//components wrapper class for EMG and UNIFORM objects

component::component(){//empty constructor
	//*****************************************************
	//Priors on Simulation for EM SEED
	//============================
	//for sigma
	alpha_0 	= 50.46;
	beta_0 		= 50.6;
	//============================
	//for lambda
	alpha_1 	= 1.823;
	beta_1 		= 0.5;
	//==============================
	//for initial length of Uniforms
	alpha_2 	= 1.297;
	beta_2 		= 19260;
	//*****************************************************
	//Priors on parameters for MAP Estimate
	ALPHA_0 = 1, BETA_0 =1; //for sigma
	ALPHA_1 = 1, BETA_1 =1; //for lambda
	ALPHA_2 = 1; //for weights, dirchlet
	ALPHA_3 = 1; //for strand probs
} 

void component::initialize(double mu, double minX, double maxX, int K, double scale, double noise_w, double noise_pi){//random seeds...
	
	if (noise_w>0){
		noise 	= NOISE(minX, maxX, noise_w, noise_pi);
		type 	= 0; 
	}else{


		//====================================
		random_device rd;
		mt19937 mt(rd());
		
		double sigma,lambda, pi_EMG, w_EMG  ;	
		double b_forward,  w_forward;
		double a_reverse,  w_reverse;

		//====================================
		//start sampling
		//for the bidirectional/EMG component
		gamma_distribution<double> dist_sigma(alpha_0,beta_0);
		gamma_distribution<double> dist_lambda(alpha_1,beta_1);
		gamma_distribution<double> dist_lengths(alpha_2,beta_2);

		sigma 		= dist_sigma(mt)/scale;
		lambda 		= dist_lambda(mt)/scale;
		b_forward 	= min(mu + (dist_lengths(mt)/scale), maxX);
		a_reverse 	= max(mu - (dist_lengths(mt)/scale),0.);

		bidir 		= EMG(mu, sigma, lambda, 1.0 / (3*K), 0.5);
		forward 	= UNI(mu, b_forward, 1.0 / (3*K), 1);
		reverse 	= UNI(a_reverse, mu, 1.0 / (3*K), -1);
		type 		= 1;
	}

} 

void component::print(){
	if (type==1){
		string text 	= bidir.print()+ "\n";
		text+=forward.print()+ "\n";
		text+=reverse.print() + "\n";

		printf(text.c_str());
	}
}

string component::write_out(){
	if (type==1){
		string text 	= bidir.print()+ "\n";
		text+=forward.print()+ "\n";
		text+=reverse.print() + "\n";

		return text;
	}
}



double component::evaluate(double x, int st){
	if (type ==0){ //this is the uniform noise component
		return noise.pdf(x, st);
	}
	if (st==1){
		bidir.ri_forward 	= bidir.pdf(x, st);
		forward.ri_forward 	= forward.pdf(x, st);
		return bidir.ri_forward + forward.ri_forward;
	}
	bidir.ri_reverse 	= bidir.pdf(x, st);
	reverse.ri_reverse 	= reverse.pdf(x, st);
		
	return bidir.ri_reverse + reverse.ri_reverse;

}


void component::add_stats(double x, double y, int st, double normalize){
	if (type==0){//noise component
		if (st==1){
			noise.r_forward+=(noise.ri_forward/normalize);
			noise.ri_forward=0;
		}else{
			noise.r_reverse+=(noise.ri_reverse/normalize);
			noise.ri_reverse=0;
		}

	}else{
		double vl, vl2;
		if (st==1){
			vl 	= bidir.ri_forward / normalize;
			vl2 = forward.ri_forward/normalize;
			bidir.ri_forward=0, forward.ri_forward=0;
			bidir.r_forward+=(vl*y);
			forward.r_forward+=(vl2*y);
		
		}else{
			vl 	= bidir.ri_reverse / normalize;
			vl2 = reverse.ri_reverse / normalize;
			bidir.ri_reverse=0, reverse.ri_reverse=0;
			bidir.r_reverse+=(vl*y);
			reverse.r_reverse+=(vl2*y);
		}
		//now adding all the conditional expections for the convolution
		double current_EY 	= bidir.EY(x, st);
		double current_EY2 	= bidir.EY2(x, st);
		double current_EX 	= x-(st*current_EY);
		bidir.ey+=current_EY*vl*y;
		bidir.ex+=current_EX*vl*y;
		bidir.ex2+=(pow(current_EX,2) + current_EY2 - pow(current_EY,2))*vl*y;	
	}
	
}

void component::reset(){
	if (type){
		bidir.ey=0, bidir.ex=0, bidir.ex2=0, bidir.r_reverse=0, bidir.r_forward=0;
		forward.r_forward=0, forward.r_reverse=0, reverse.r_reverse=0, reverse.r_forward=0;
	}else{
		noise.r_forward=0,noise.r_reverse=0;
	}
}

double component::get_all_repo(){
	if (type==1){
		return bidir.r_forward+bidir.r_reverse+forward.r_forward+reverse.r_reverse;
	}
	return noise.r_forward+noise.r_reverse;

}

void component::update_parameters(double N, int K){
	if (type==1){
		//first for the bidirectional
		double r 	= bidir.r_forward + bidir.r_reverse;
		bidir.pi 	= (bidir.r_forward + ALPHA_3) / (r + ALPHA_3*2);
		bidir.w 	= (r + ALPHA_2) / (N + ALPHA_2*K*3 + K*3) ;
		bidir.mu 	= bidir.ex / (r+0.001);

		
		bidir.si 	= pow(abs((1. /(r + 3 + ALPHA_0 ))*(bidir.ex2-2*bidir.mu*bidir.ex + r*pow(bidir.mu,2) + 2*BETA_0  )), 0.5);
		bidir.l 	= min((r+ALPHA_1) / (bidir.ey + ALPHA_1), 2.);
		//now for the forward and reverse strand elongation components
		forward.w 	= forward.r_forward / N;
		reverse.w 	= reverse.r_reverse / N;
		forward.a 	= bidir.mu, reverse.b=bidir.mu;

	}
}

bool component::check_elongation_support(){
	if (forward.b <=forward.a and bidir.mu==0){
		return true;
	}
	else if(reverse.b <= reverse.a and bidir.mu==0){
		return true;
	}
}





//=========================================================
//FIT function this where we are running EM rounds
//=========================================================
//helper functions for fit

double sum(double * X, int N){
	double vl=0;
	for (int i = 0; i < N; i++){
		vl=(vl+X[i]);
	}
	return vl;

}
double LOG(double x){
	if (x <= 0){
		return nINF;
	}
	return log(x);
}


double calc_log_likelihood(component * components, int K, segment * data){
	double ll 	= 0;
	double forward, reverse;
	for (int i = 0 ; i < data->XN; i++){
		forward=0, reverse=0;
		for (int k = 0; k < K; k++){
			forward+=(components[k].evaluate(data->X[0][i], 1));
			reverse+=(components[k].evaluate(data->X[0][i], -1));
		}
		ll+=LOG(forward)*data->X[1][i] + LOG(reverse)*data->X[2][i];
	}
	return ll;
}	

double move_uniforom_support(component * components, int K, int add, segment * data, double move, double base_ll){
	//===========================================================================
	//normal distribution centered around 0 and some variance, how much to move 
	//uniform supports
	random_device rd;
	mt19937 mt(rd());
	normal_distribution<double> dist_uni(0,move);
	vector<int> components_that_made_it;
	vector<int> type_that_made_it;
	vector<double> and_their_moves;
	double ll, step;
	double maxLL 	= base_ll;
	for (int k = 0; k < K;k++){
		//try the forward
		step 	= dist_uni(mt) ;
		components[k].forward.delta_b=step;
		ll 	 	= calc_log_likelihood(components, K+add, data);
		if (ll > base_ll){
			components_that_made_it.push_back(k);
			type_that_made_it.push_back(1);
			and_their_moves.push_back(step);
		}
		components[k].forward.delta_b=0;
		step 	= dist_uni(mt);;
		components[k].reverse.delta_a=step ;
		ll 		= calc_log_likelihood(components, K+add, data);
		if (ll > base_ll){
			components_that_made_it.push_back(k);
			type_that_made_it.push_back(-1);
			and_their_moves.push_back(step);
		}		
		components[k].reverse.delta_a=0;
	}
	for (int u = 0 ; u < components_that_made_it.size(); u++){
		if (type_that_made_it[u]==1){
			components[components_that_made_it[u]].forward.b+=and_their_moves[u];
		}else{
			components[components_that_made_it[u]].reverse.a+=and_their_moves[u];	
		}
	}
	ll 		= calc_log_likelihood(components, K+add, data);
	return ll;
}

//=========================================================
//For Classifier class / wrapper around EM

classifier::classifier(int k, double ct, int mi, double nm,
	double move, double R_MU){
	K 						= k ;
	seed 					= true;
	convergence_threshold 	= ct;
	max_iterations 			= mi;
	noise_max 				= nm;
	p 						= 0.8;
	move 					= move;
	resets 					= 0;
	last_diff 				= 0;
	r_mu 					= R_MU;

}

classifier::classifier(){};//empty constructor

int classifier::fit(segment * data, vector<double> mu_seeds){
	//=========================================================================
	//for resets
	random_device rd;
	mt19937 mt(rd());
	uniform_real_distribution<double> dist_uni(data->minX,data->maxX);
	

	int i;
	double l 	= data->maxX - data->minX;
	pi 	= sum(data->X[1], data->XN)/ data->N;
	double vl 	= 1.0 / l;
	if (K==0){
		//calc_likeihood coming from uniform model, only
		ll 	= 0;
		for (int i = 0; i < data->XN; i ++){
			ll+=(LOG(vl*(pi) )*data->X[1][i]);
			ll+=(LOG(vl*(1-pi))*data->X[2][i]);
		}

		return 1;
	}
	int add 	= noise_max>0;
	components 	= new component[K+add];
	//===========================================================================
	//random seeds, initialize
	double mu;
	for (int k = 0; k < K; k++){

		if (mu_seeds.size()>0){
			i 	= sample_centers(mu_seeds ,  p);
			mu 	= mu_seeds[i];
			if (r_mu > 0){
				normal_distribution<double> dist_r_mu(mu, r_mu);
				mu 		= dist_r_mu(mt);
			}
		}else{
			mu 			= dist_uni(mt);
		}
		components[k].initialize(mu, data->minX, data->maxX, K, data->SCALE , 0., 0.);
		mu_seeds.erase (mu_seeds.begin()+i);	
	}
	if (add){
		components[K].initialize(0., data->minX, data->maxX, 0., 0. , noise_max, pi);
	}
		
		
	//===========================================================================
	int t 			= 0; //EM loop ticker
	double prevll 	= nINF; //previous iterations log likelihood
	converged 		= false; //has the EM converged?
	double norm_forward, norm_reverse,N; //helper variables
	while (t < max_iterations && not converged){
		
		//******
		//reset old sufficient statistics
		for (int k=0; k < K+add; k++){
			components[k].reset();
		}
		//******
		//E-step
		for (int i =0; i < data->XN;i++){
			norm_forward=0;
			norm_reverse=0;
			for (int k=0; k < K+add; k++){
				if (data->X[1][i]){//if there is actually data point here...
					norm_forward+=components[k].evaluate(data->X[0][i],1);
				}
				if (data->X[2][i]){//if there is actually data point here...
					norm_reverse+=components[k].evaluate(data->X[0][i],-1);
				}
			}
			//now we need to add the sufficient statistics
			for (int k=0; k < K+add; k++){
				if (norm_forward){
					components[k].add_stats(data->X[0][i], data->X[1][i], 1, norm_forward);
				}
				if (norm_reverse){
					components[k].add_stats(data->X[0][i], data->X[2][i], -1, norm_reverse);
				}
			}

		
		}
		//******
		//M-step
		//get normalizing constant
		N=0;
		for (int k = 0; k < K+add; k++){
			N+=(components[k].get_all_repo());
		}
		//update the new parameters
		for (int k = 0; k < K+add; k++){
			components[k].update_parameters(N, K);
		}
		

		ll 	= calc_log_likelihood(components, K+add, data);
		//******
		//Move Uniform support		
		ll 	= move_uniforom_support(components, K, add, data,move, ll);
		if (abs(ll-prevll)<convergence_threshold){
			converged=true;
		}
		last_diff=abs(ll-prevll);
		prevll=ll;

		t++;
	}
	return 1;
}
string classifier::print_out_components(){
	string text 	= "";
	for (int k = 0; k < K; k++){
		text+=components[k].write_out();
	}
	return text;
}














