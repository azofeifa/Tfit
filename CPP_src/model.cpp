#include <mpi.h>
#include "model.h"
#include "load.h"
#include "template_matching.h"
#include <math.h> 
#include <limits>
#include <iostream>
#include <algorithm>
#include <unistd.h>
#include <random>
#include "omp.h"


//=============================================
//Helper functions
double IN(double x){ //Standard Normal PDF 
	return exp(-pow(x,2)*0.5)/sqrt(2*M_PI);
	
}
double IC(double x){ //Standard Normal CDF
	return 0.5*(1+erf(x/sqrt(2)));
}


double R(double x){ //Mills Ratio
	if (x > 4){
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
		return (w*pi) / abs(b-a);
	}
	return (w*(1-pi)) / abs(b-a);
}
NOISE::NOISE(){}
NOISE::NOISE(double A, double B, double W, double PI){
	a=A;
	b=B;
	w=W;
	pi=PI;
}

//=============================================
//Normal Distribution Class

NORMAL::NORMAL(){}
NORMAL::NORMAL(double MU, double SI, double PI, double W){
	mu=MU,si=SI,pi=PI,w=W;
}
double NORMAL::pdf(double x, int s){
	if ((s == 1 and pi==1) or (s== -1 and pi == 0)){
		return (w/(si*sqrt(2*M_PI))  )*exp(-pow(x-mu,2)/(2*pow(si,2)));
	}else{
		return 0;
	}

}

//=============================================
//Uniform Class
double UNI::pdf(double x, int strand){ //probability density function
	double p;
	if (w==0){
		return 0;
	}

	if ( a<= x and x <=b){

		p= w / abs(b- a);
		p= p*pow(pi, max(0, strand) )*pow(1.-pi, max(0, -strand) );
		return p;
	}
	return 0;
}

string UNI::print(){
	string text = ("U: " + to_string(a) + "," + to_string(b) 
	+ "," + to_string(w) + "," + to_string(pi));
	return text;
}	
UNI::UNI(double start, double stop, double w_i, int strand, int POS, double Pi){
	a 		= start;
	b 		= stop;
	w 		= w_i;
	st 		= strand;
	pos 	= POS;
	if (st==1){
		pi=1;
	}else{
		pi=0;
	}
	//===================
	//this oversets the constraint that uniform must take either 
	//forward or reverse data points
	pi 		= Pi;
	//===================
	delta_a=0;
	delta_b=0;
	ri_forward=0, ri_reverse=0;


}

UNI::UNI(){} //empty constructor


//=============================================
//Exponentially Modified Gaussian class

string EMG::print(){
	string text 	= ("N: " + to_string(mu)+","+to_string(si)
		+ "," + to_string(l) + "," + to_string(w) + "," + to_string(pi) + "," + to_string(foot_print) );
	return text;
}

EMG::EMG(){
	foot_print=0;
}//empty constructor

EMG::EMG(double MU, double SI, double L, double W, double PI ){
	mu 	= MU;
	si 	= SI;
	l  	= L;
	w 	= W;
	pi 	= PI;
}


double EMG::pdf(double z, int s ){
	if (w==0){
		return 0;
	}
	if (s==1){
		z-=foot_print;
	}else{
		z+=foot_print;
	}
	double vl 		= (l/2)*(s*2*(mu-z) + l*pow(si,2));
	double p;
	if (vl > 100){ //potential for overflow, inaccuracies
		p 			= l*IN((z-mu)/si)*R(l*si - s*((z-mu)/si));
	}else{
		p 			= (l/2)*exp(vl)*erfc((s*(mu-z) + l*pow(si ,2) )/(sqrt(2)*si));
	}
	vl 				= p*w*pow(pi, max(0, s) )*pow(1.-pi, max(0, -s) );
	if (checkNumber(vl)){
		return vl;
	}
	return 0.;

}
double EMG::EY(double z, int s){
	if (s==1){
		z-=foot_print;
	}else{
		z+=foot_print;
	}
	
	return max(0. , s*(z-mu) - l*pow(si, 2) + (si / R(l*si - s*((z-mu)/si))));
}
double EMG::EY2(double z, int s){
	if (s==1){
		z-=foot_print;
	}else{
		z+=foot_print;
	}
	
	return pow(l,2)*pow(si,4) + pow(si, 2)*(2*l*s*(mu-z)+1 ) + pow(mu-z,2) - ((si*(l*pow(si,2) + s*(mu-z)))/R(l*si - s*((z-mu)/si) ));
	
}

//=========================================================
//components wrapper class for EMG and UNIFORM objects

component::component(){//empty constructor
	foot_print 			= 0;
	forward_neighbor 	= NULL;
	reverse_neighbor 	= NULL;
} 

void component::set_priors(double s_0, double s_1, 
	double l_0, double l_1, double w_0,double strand_0, double N, int K){
	//============================
	//for sigma
	alpha_0 	= 20.46;
	beta_0 		= 10.6;
	//============================
	//for lambda
	alpha_1 	= 20.823;
	beta_1 		= 0.5;
	//==============================
	//for initial length of Uniforms
	alpha_2 	= 1.297;
	beta_2 		= 8260;

	//*****************************************************
	//Priors on parameters for MAP Estimate
	ALPHA_0 = s_0, BETA_0 =s_1; //for sigma
	ALPHA_1 = l_0, BETA_1 =l_1; //for lambda
	ALPHA_2 = w_0; //for weights, dirchlet
	ALPHA_3 = strand_0; //for strand probs
	//bidir.w 	= (r + ALPHA_2) / (N + ALPHA_2*K*3 + K*3) ;
		
	w_thresh= ( ALPHA_2 ) / (N + ALPHA_2*K*3 + K*3 );
	
}

bool check_uniform_support(component c, int forward){
	if (forward==1){
		if (c.forward.b < (c.bidir.mu + c.bidir.si + (1.0 / c.bidir.l) )){
			return false;
		}
		return true;
	}
	else{
		if (c.reverse.a < (c.bidir.mu - c.bidir.si - (1.0 / c.bidir.l) )){
			return false;
		}
		return true;

	}
}
int get_nearest_position(segment * data, double center, double dist){
	int i;

	if (dist < 0 ){
		i=0;
		while (i < (data->XN-1) and (data->X[0][i] -center) < dist){
			i++;
		}
	}else{
		i=data->XN-1;
		while (i >0 and (data->X[0][i] - center) > dist){
			i--;
		}
	}
	return i;
}


void component::initialize(double mu, segment * data , int K, double scale, double noise_w, 
	double noise_pi, double fp){//random seeds...
	foot_print 	= fp;
	EXIT=false;
	if (noise_w>0){
		noise 	= NOISE(data->minX, data->maxX, 
			noise_w, 0.5);
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
		uniform_real_distribution<double> dist_lambda_2(1, 500);
		uniform_real_distribution<double> dist_sigma_2(1, 50);
		uniform_real_distribution<double> dist_footprint(0, 100);
		
		gamma_distribution<double> dist_lengths(1,( (data->maxX-data->minX)/(K)));
		
		sigma 		= dist_sigma_2(mt)/scale;
		lambda 		= scale/dist_lambda_2(mt) ;
		double dist = data->maxX - mu+(1.0/lambda) ;
		int j 		= get_nearest_position(data, mu, dist);
		
		forward 	= UNI(mu+(1.0/lambda), data->maxX, 1.0 / (3*K), 1, j, 0.5);
		
		dist 		= data->minX - mu - (-1.0/lambda) ;
		j 			= get_nearest_position(  data, mu, dist);
		
		a_reverse 			= data->X[0][j];

		bidir 				= EMG(mu, sigma, lambda, 1.0 / (3*K), 0.5);
		bidir.foot_print 	= 50/100.;
		bidir.prev_mu 		= mu;
		reverse 	= UNI(data->minX, mu-(1.0/lambda), 1.0 / (3*K), -1, j,0.5);
		type 		= 1;
		
	}

} 

void component::initialize_bounds(double mu, segment * data , int K, double scale, double noise_w, 
	double termination, double fp, double forward_bound, double reverse_bound){//random seeds...
	foot_print 	= fp;
	EXIT=false;
	if (noise_w>0){
		noise 	= NOISE(data->minX, data->maxX, 
			noise_w, 0.5);
		type 	= 0; 
	}else{

		int complexity=1;
		if (data->strand == "."){
			complexity=3;
		}else{
			complexity=2;
		}

		//====================================
		random_device rd;
		mt19937 mt(rd());
		
		double sigma,lambda, pi_EMG, w_EMG  ;	
		double b_forward,  w_forward;
		double a_reverse,  w_reverse;

		//====================================
		//start sampling
		//for the bidirectional/EMG component
		uniform_real_distribution<double> dist_lambda_2(1, 250);
		uniform_real_distribution<double> dist_sigma_2(1, 250);
		gamma_distribution<double> dist_lengths(1,( (data->maxX-data->minX)/(K)));
		
		sigma 				= dist_sigma_2(mt)/scale;
		lambda 				= scale/dist_lambda_2(mt) ;
		double dist 		=  (1.0/lambda);
		int j 				= get_nearest_position(data, mu, dist);
		int k 				= get_nearest_position(data, mu, forward_bound-mu);
		double pi 			= 0.5;
		if (data->strand == "+"){
			pi 				= 1.0;
		}else if (data->strand=="-") {
			pi 				= 0.;
		}
		forward 			= UNI(mu+(1.0/lambda), data->maxX, 1.0 / (complexity*K), 1, j, pi);
		forward.j=j, forward.k=k;
		

		bidir 				= EMG(mu, sigma, lambda, 1.0 / (complexity*K), 0.5);
		bidir.foot_print 	= fp;
		dist 				=  -(1.0/lambda);
		j 					= get_nearest_position(  data, mu, dist);
		k 					= get_nearest_position(data, mu, reverse_bound-mu);
		reverse 			= UNI(data->minX, mu-(1.0/lambda),  1.0 / (complexity*K) , -1, j,1-pi);
		reverse.j=k, reverse.k=j;
		termination 		= (termination > 1);
		
		type 		= 1;
		
	}

} 



void component::print(){
	if (type==1){
		string text 	= bidir.print()+ "\n";
		text+=forward.print()+ "\n";
		text+=reverse.print() + "\n";
		cout<<text;
	}else{
		cout<<"NOISE: " << noise.w<<"," <<noise.pi<<endl;
	}
}

string component::write_out(){
	if (type==1){
		string text 	= bidir.print()+ "\n";
		text+=forward.print()+ "\n";
		text+=reverse.print() + "\n";

		return text;
	}
	return "";
}


void sort_components(component components[], int K){
	bool sorted=true;
	while (sorted){
		sorted=false;
		for (int k = 0; k < K-1; k++){
			if (components[k].bidir.mu > components[k+1].bidir.mu){
				component cp 		= components[k];
				components[k] 		= components[k+1];
				components[k+1] 	= cp;
				sorted 				= true;
			}
		}
	}
}


double component::evaluate(double x, int st){
	if (type ==0){ //this is the uniform noise component
		return noise.pdf(x, st);
	}
	if (st==1){
		bidir.ri_forward 	= bidir.pdf(x, st);
		forward.ri_forward 	= forward.pdf(x, st);
		reverse.ri_forward 	= reverse.pdf(x,st);
		return bidir.ri_forward + forward.ri_forward + reverse.ri_forward;
	}
	bidir.ri_reverse 	= bidir.pdf(x, st);
	reverse.ri_reverse 	= reverse.pdf(x, st);
	forward.ri_reverse 	= forward.pdf(x, st);
	return bidir.ri_reverse + reverse.ri_reverse + forward.ri_reverse;
}
double component::pdf(double x, int st){
	if (type==0){
		return noise.pdf(x,st);
	}
	if (st==1){
		return bidir.pdf(x,st) + forward.pdf(x, st);
	}
	return bidir.pdf(x,st) + reverse.pdf(x,st);
}


void component::add_stats(double x, double y, int st, double normalize){
	if (type==0){//noise component
		if (st==1){
			noise.r_forward+=(y*noise.ri_forward/normalize);
			noise.ri_forward=0;
		}else{
			noise.r_reverse+=(y*noise.ri_reverse/normalize);
			noise.ri_reverse=0;
		}

	}else{
		double vl, vl2, vl3;
		if (st==1){
			vl 	= bidir.ri_forward / normalize;
			vl2 = forward.ri_forward/normalize;
			vl3 = reverse.ri_forward/normalize;
			bidir.ri_forward=0, forward.ri_forward=0;
			bidir.r_forward+=(vl*y);
			forward.r_forward+=(vl2*y);
			reverse.r_forward+=(vl3*y);
		
		}else{
			vl 	= bidir.ri_reverse / normalize;
			vl2 = reverse.ri_reverse / normalize;
			vl3 = forward.ri_reverse / normalize;
			bidir.ri_reverse=0, reverse.ri_reverse=0;
			bidir.r_reverse+=(vl*y);

			reverse.r_reverse+=(vl2*y);
			forward.r_reverse+=(vl3*y);
		}
		//now adding all the conditional expections for the convolution
		if (vl > 0 and y > 0){
			double current_EY 	= bidir.EY(x, st);
			double current_EY2 	= bidir.EY2(x, st);
			double current_EX 	= x-(st*current_EY)-bidir.foot_print*st;
			//	self.C+=max( ((z-self.mu) -E_Y) *r,0)
			// 	self.C+=max((-(z-self.mu) -E_Y)   *r ,0)
			bidir.C+=max((st*(x-bidir.mu) - current_EY  )*vl*y,0.0);
			bidir.ey+=current_EY*vl*y;
			bidir.ex+=current_EX*vl*y;
			bidir.ex2+=(pow(current_EX,2) + current_EY2 - pow(current_EY,2))*vl*y;	
		}
	}
	
}

void component::reset(){
	if (type){
		bidir.C=0;
		bidir.ey=0, bidir.ex=0, bidir.ex2=0, bidir.r_reverse=0, bidir.r_forward=0;
		bidir.ri_forward=0, forward.ri_forward=0, forward.ri_reverse=0;
		bidir.ri_reverse=0, reverse.ri_reverse=0, reverse.ri_forward=0;
		forward.r_forward=0, forward.r_reverse=0, reverse.r_reverse=0, reverse.r_forward=0;
		forward.delta_a=0, forward.delta_b=0, reverse.delta_a=0, reverse.delta_b=0;
	}else{
		noise.r_forward=0,noise.r_reverse=0;
		noise.ri_reverse=0,noise.ri_forward=0 ;
		
	}
}

double component::get_all_repo(){
	if (type==1){
		return bidir.r_forward+bidir.r_reverse+forward.r_forward+reverse.r_reverse;
	}
//	return noise.r_forward+noise.r_reverse;
	return 0.;

}
double get_sum(segment * data, int j, int k, int st){
	double S 	= 0;
	for (int i = j; i <k;i++){
		S+=data->X[st][i];
	}
	return S;
}


void update_j_k( component * components,segment * data, int K, double N){
	for (int k = 0 ; k < K;k++){
		if (k > 0){
			components[k].reverse_neighbor 	= &components[k-1]; 
		}
		if (k+1 < K){
			components[k].forward_neighbor 	= &components[k+1];
		}
	}

	double center, w_thresh;

	for (int k = 0; k < K; k++){
		int j 	= k;
		//for the forward
		w_thresh= ( components[k].ALPHA_2 ) / (N + components[k].ALPHA_2*K*3 + K*3 );
		while ( j < K and components[j].forward_neighbor!=NULL and  components[j].forward_neighbor->bidir.w  < pow(10,-2)   ){
			j++;
		}
		double delta 				= (1.0 / components[k].bidir.l);
		components[k].forward.j 	= get_nearest_position(data, components[k].bidir.mu, delta  );
		
		if (j < K and components[j].forward_neighbor!=NULL  ){
			center 						= components[k].bidir.mu + delta;
			components[k].forward.k 	= get_nearest_position(data, center, components[j].forward_neighbor->bidir.mu -center  );
		}else{
			components[k].forward.k 	= get_nearest_position(data, components[k].bidir.mu,data->maxX -components[k].bidir.mu );			
		}
		if (components[k].forward.j >= components[k].forward.k ){
			components[k].forward.j 	= components[k].forward.k;
			components[k].forward.w 	= 0;
		}
		
		j 	= k;
		while (j >=0 and components[j].reverse_neighbor!=NULL and  components[j].reverse_neighbor->bidir.w   < pow(10,-2) ){
			j--;
		}



		if (K>j and j>=0 and components[j].reverse_neighbor!=NULL){
			center 						= components[k].bidir.mu-delta;
			components[k].reverse.j 	= get_nearest_position(data, center, components[j].reverse_neighbor->bidir.mu - center );
		}else{
			components[k].reverse.j 	= get_nearest_position(data, components[k].bidir.mu, data->minX-components[k].bidir.mu );			
		}
		components[k].reverse.k 		= get_nearest_position(data, components[k].bidir.mu, -delta);
		if (components[k].reverse.j >= components[k].reverse.k ){
			components[k].reverse.j 	= components[k].reverse.k;
			components[k].reverse.w 	= 0;
		}		
	}
}
void update_l(component * components, segment * data, int K){
	for (int k 	= 0; k < K; k++){
		//forward
		double left_SUM=0, right_SUM=get_sum(data,components[k].forward.j ,components[k].forward.k,1);
		double N 		= left_SUM+right_SUM;
		double null_vl 	= 0.05 / (data->maxX - data->minX  );
		double null_ll 	= N*LOG(1. / (data->maxX - data->minX  ));
		double null_BIC = -2*null_ll + LOG(N);
		double vl 		= 0, w=0;
		double mod_ll, mod_BIC;
		double prev_prev=0, prev=0, current=0;
		double BIC_best = 0;
		int arg_l 	= components[k].forward.k;

		for (int l = components[k].forward.j; l < components[k].forward.k; l++ ){
			left_SUM+=data->X[1][l];
			right_SUM-=data->X[1][l];
			vl 		= 1.0/(data->X[0][l]-data->X[0][components[k].forward.j]);
			w 		= left_SUM/(N);
			mod_ll 	= LOG(w*vl)*left_SUM + LOG(null_vl)*right_SUM ;
			mod_BIC = -2*mod_ll + 5*LOG(N);
			current = null_BIC/mod_BIC;
			if (prev > current and prev > prev_prev and prev > 1.0){
				if (prev > BIC_best){
					BIC_best=prev, arg_l=l;
				}
			}
			prev_prev=prev;
			prev 	= current;			
		}
		components[k].forward.b 	= data->X[0][arg_l];
		//reverse
		arg_l 	= components[k].reverse.j;
		left_SUM=0, right_SUM=get_sum(data,components[k].reverse.j,components[k].reverse.k,2 );
		N 		= left_SUM+right_SUM;
		null_ll 	= N*LOG(1. / (data->maxX - data->minX  ));
		null_BIC = -2*null_ll + LOG(N);
		prev_prev=0, prev=0, current=0,BIC_best = 0;
		for (int l = components[k].reverse.j; l < components[k].reverse.k; l++ ){
			left_SUM+=data->X[2][l];
			right_SUM-=data->X[2][l];
			vl 		= 1.0/(data->X[0][components[k].reverse.k] - data->X[0][l]);
			w 		= right_SUM/(N);
			mod_ll 	= LOG(null_vl)*left_SUM + LOG(w*vl)*right_SUM ;
			mod_BIC = -2*mod_ll + 5*LOG(N);
			current = null_BIC/mod_BIC;
			if (prev > current and prev > prev_prev and prev > 1.0){
				if (prev > BIC_best){	
					BIC_best=prev, arg_l=l;
				}
			}
			prev_prev=prev;
			prev 	= current;
		}
		components[k].reverse.a 	= data->X[0][arg_l];
	}
}



void component::update_parameters(double N, int K){
	if (type==1){
		//first for the bidirectional
		double r 	= bidir.r_forward + bidir.r_reverse;
		bidir.pi 	= (bidir.r_forward + ALPHA_3) / (r + ALPHA_3*2);
		bidir.w 	= (r + ALPHA_2) / (N + ALPHA_2*K*3 + K*3) ;
		bidir.mu 	= bidir.ex / (r+0.001);

		
		bidir.si 	= pow(abs((1. /(r + 3 + ALPHA_0 ))*(bidir.ex2-2*bidir.mu*bidir.ex + 
			r*pow(bidir.mu,2) + 2*BETA_0  )), 0.5) ;
		if (bidir.si > 50){
			EXIT 	= true;
			bidir.w = 0;
		}
		if ((r / N) < pow(10,-5) ){
			EXIT 	= true;
		}
		bidir.l 	= min((r+ALPHA_1) / (bidir.ey + BETA_1), 5.);
		bidir.l 	= max(bidir.l, 0.05);
		if (bidir.l < 0.05  ){
			EXIT 	= true;
			bidir.w = 0;
		}
		if (abs(bidir.mu-bidir.prev_mu)< 0.001 ){
			bidir.move_fp 	= true;
		}
		else{
			bidir.prev_mu 	= bidir.mu;
		}
		if (bidir.move_fp){

			bidir.foot_print 	= min( max(bidir.C / (r+0.1),0.0) , 5.0);
			bidir.foot_print 	= floor((bidir.foot_print*10000.0))/10000.0;
		}
		//bidir.foot_print 	= 0.0;
		//now for the forward and reverse strand elongation components
		forward.w 	= (forward.r_forward + ALPHA_2) / (N+ ALPHA_2*K*3 + K*3);
		reverse.w 	= (reverse.r_reverse + ALPHA_2) / (N+ ALPHA_2*K*3 + K*3);
		forward.a 	= bidir.mu  , reverse.b=bidir.mu ;

		//update PIS, this is obviously overwritten if we start the EM seeder with 0/1
		forward.pi 	= (forward.r_forward + 1) / (forward.r_forward + forward.r_reverse+2);
		reverse.pi 	= (reverse.r_forward + 1)/ (reverse.r_forward + reverse.r_reverse+2);
		if (bidir.w==0){
			forward.w 	= 0;
			reverse.w 	= 0;
		}
	}
}

void check_mu_positions(component * components, int K){
	double dist 	= 0;
	double dist_a  	= 0;
	double dist_b 	= 0;
	double dist_c 	= 0;
	for (int i = 0; i < K; i++){
		for (int j = i+1; j < K; j++){
			dist_a 	= components[i].bidir.mu + components[i].bidir.si + (1.0/components[i].bidir.l);
			dist_b 	= components[j].bidir.mu - components[j].bidir.si - (1.0/components[j].bidir.l);
			if (dist_b < dist_a   ) {
				components[i].EXIT 	= true, components[j].EXIT=true;
			}
		}
	}
}


bool component::check_elongation_support(){
	if (forward.b <=forward.a and bidir.mu==0){
		return true;
	}
	else if(reverse.b <= reverse.a and bidir.mu==0){
		return true;
	}
	return false;
}


void component::initialize_with_parameters(vector<double> init_parameters, segment * data, 
	int K, double left, double right, double fp){
	foot_print 	= fp;
	if (init_parameters.empty()){
		noise 	= NOISE(data->minX, data->maxX, 0.01, 0.5);
		type 	= 0; 
	}else{
		double mu 	= init_parameters[0];
		double si 	= init_parameters[1];
		double l 	= init_parameters[2];
		double pi 	= init_parameters[3];
		bidir 		= EMG(mu, si, l, 1.0 / (3*K), pi);//bidir component
		//now choose supports of forward and reverse
		random_device rd;
		mt19937 mt(rd());
		uniform_real_distribution<double> left_move( left , mu + (1.0/l) );
		uniform_real_distribution<double> right_move( mu-(1.0/l ), right );
			
		forward 	= UNI(mu+(1.0/l), right_move(mt) , 1. / double(K), 1, 0, 1);
				
		reverse 	= UNI(left_move(mt), mu-(1.0/l ), 1. / double(K), -1, 0,0.);
		type 		= 1;
	}	
}

void component::initialize_with_parameters2(vector<double> init_parameters, segment * data, 
	int K, double left, double right,double fp){
	foot_print 	= fp;
	if (init_parameters.empty()){
		noise 	= NOISE(data->minX, data->maxX, 0.01, 0.5);
		type 	= 0; 
	}else{
		double mu 	= init_parameters[0];
		double si 	= init_parameters[1];
		double l 	= init_parameters[2];
		double pi 	= init_parameters[3];
		bidir 		= EMG(mu, si, l, 1.0 / (3*K), pi);//bidir component
		//now choose supports of forward and reverse
		
		forward 	= UNI(mu+(1.0/l), right , 1. / double(3*K), 1, 0, 1);
				
		reverse 	= UNI(left, mu-(1.0/l ), 1. / double(3*K), -1, 0,0.);
		type 		= 1;
	}	
}


//=========================================================
//FIT function this where we are running EM rounds
//=========================================================
//helper functions for fit

double sum(double * X, int N){
	double vl=0;
	for (int i = 0; i < N; i++){
		vl+=X[i];
	}
	return vl;

}
double LOG(double x){
	if (x <= 0){
		return nINF;
	}
	return log(x);
}
double pdf_all(component * components, double x, double f_y, double r_y,int K){
	double forward=0;
	double reverse=0;
	for (int k = 0; k < K; k++){
		forward+=(components[k].pdf(x, 1));
		reverse+=(components[k].pdf(x, -1));
	}
	return 	(LOG(forward)*f_y + LOG(reverse)*r_y); 

}

double calc_log_likelihood(component * components, int K, segment * data){
	double ll 	= 0;
	int N 		= data->XN;
	#pragma omp parallel for reduction(+:ll)
	for (int i = 0 ; i < N; i++){
		ll+=pdf_all(components, data->X[0][i], data->X[1][i], data->X[2][i], K) ; 
	}
	return ll;
}	
int get_direction(uniform_int_distribution<int> direction,mt19937 mt){
	if (direction(mt)==0){
		return 1;
	}
	return -1;
}
int get_new_position(geometric_distribution<int> dist_uni, mt19937 mt, 
	int pos, int N, int direction, int st, segment * data){
	int ct 	= dist_uni(mt)+1;
	int i = pos;
	int j = 0;

	if (direction==1){
		i 	= min(N-1, i+ct);
	}else{
		i 	= max(0, i-ct);
	}
	//printf("old: %d,new: %d, change: %d,old: %f, new: %f\n", pos,i, ct*direction, data->X[0][pos], data->X[0][i] );
	return i;
}

double move_uniforom_support(component * components, int K, int add, 
	segment * data, double move, double base_ll){
	//===========================================================================
	//normal distribution centered around 0 and some variance, how much to move 
	//uniform supports
	random_device rd;
	mt19937 mt(rd());
	geometric_distribution<int> dist_uni(0.9);
	uniform_int_distribution<int> direction(0,1);
	int 	steps[K][2];
	double  new_bounds[K][2];
	double ll;
	double prev_a, prev_b;
	for (int k = 0; k < K; k++){
		prev_b=components[k].forward.b, prev_a=components[k].reverse.a;
		steps[k][0] 	= get_new_position(dist_uni, mt, 
			components[k].forward.pos, data->XN, get_direction(direction, mt), 1, data);
		steps[k][1] 	= get_new_position(dist_uni, mt, 
			components[k].reverse.pos, data->XN, get_direction(direction, mt), 2, data);
		//firs the forward
		components[k].forward.b=data->X[0][steps[k][0]];
		if (check_uniform_support(components[k], 1)){
			ll 	= calc_log_likelihood(components, K+add, data);
			//printf("%f,%f,%f,%f\n", prev_b, components[k].forward.b, base_ll,ll );
			if (ll > base_ll){
				new_bounds[k][0] 	= components[k].forward.b;
			}else{
				new_bounds[k][0] 	= prev_b;
				steps[k][0] 		= components[k].forward.pos;
			}
		}else{
			new_bounds[k][0] 	= prev_b;
			steps[k][0] 		= components[k].forward.pos;
		}
		components[k].forward.b 	= prev_b;
		//now reverse
		components[k].reverse.a=data->X[0][steps[k][1]];
		if (check_uniform_support(components[k], 0)){

			ll 	= calc_log_likelihood(components, K+add, data);
			if (ll > base_ll){
				new_bounds[k][1] 	= components[k].reverse.a;
			}else{
				new_bounds[k][1] 	= prev_a;
				steps[k][1] 		= components[k].reverse.pos;
			}
			//printf("%f,%f,%f,%f\n", prev_a, components[k].reverse.a, base_ll,ll );
		}else{
				new_bounds[k][1] 	= prev_a;
				steps[k][1] 		= components[k].reverse.pos;			
		}
		components[k].reverse.a 	= prev_a;
	}
	for (int k =0; k < K;k++){
		components[k].forward.pos 	= steps[k][0];
		components[k].forward.b 	= new_bounds[k][0];
		components[k].reverse.pos 	= steps[k][1];
		components[k].reverse.a 	= new_bounds[k][1];
 		
			
	}



	ll 		= calc_log_likelihood(components, K+add, data);
	return ll;
}

void update_weights_only(component * components, segment * data, int K, int add){
	
	double EX[K+add][3];
	double current[K+add][3][2];
	double x,f,r, normed_f, normed_r, N,W;
	for (int k =0; k < K+add; k++){
		EX[k][0]=0,EX[k][1]=0,EX[k][2]=0;
		current[K+add][0][0]=0, current[K+add][1][0]=0, current[K+add][2][0]=0;

	}
	for (int i =0; i< data->XN; i++){
		x 	= data->X[0][i], f 	=data->X[1][i], r=data->X[2][i];
		normed_f = 0, normed_r 	= 0;
		for (int k = 0; k<K; k++){
			current[k][0][0] 	= components[k].bidir.pdf(x,1);
			current[k][1][0] 	= components[k].forward.pdf(x,1);
			current[k][2][0] 	= components[k].reverse.pdf(x,1);
			current[k][0][1] 	= components[k].bidir.pdf(x,-1);			
			current[k][1][1] 	= components[k].forward.pdf(x,-1);
			current[k][2][1] 	= components[k].reverse.pdf(x,-1);
			
			normed_f+=(current[k][0][0]+ current[k][1][0] + current[k][2][0]);
			normed_r+=(current[k][0][1]+ current[k][1][1] + current[k][2][1]);	
		}
		if (add){
			current[K][0][0] 	= components[K].noise.pdf(x,1);
			current[K][0][1] 	= components[K].noise.pdf(x,-1);
			normed_f+=current[K][0][0];
			normed_r+=current[K][0][1];			
		}
		for (int k =0; k < K+add; k++){

			if (normed_f){
				if (k < K){
					EX[k][0]+=f*(current[k][0][0]/normed_f);
					EX[k][1]+=f*(current[k][1][0]/normed_f);
					EX[k][2]+=f*(current[k][2][0]/normed_f);
				}else{
					EX[k][0]+=f*(current[k][0][0]/normed_f);
				}
			}
			if (normed_r){
				if (k < K){
					EX[k][0]+=r*(current[k][0][1]/normed_r);
					EX[k][1]+=r*(current[k][1][1]/normed_r);
					EX[k][2]+=r*(current[k][2][1]/normed_r);
				}else{
					EX[k][0]+=r*(current[k][0][1]/normed_r);		
				}
			}
		}
	}
	N 	= 0;
	W 	= 0;
	for (int k = 0; k < K+add; k++){
		N+=(EX[k][0]+EX[k][1]+EX[k][2]);
	}
	for (int k = 0; k < K; k++){
		components[k].bidir.w 		= EX[k][0] / N;
		components[k].forward.w 	= EX[k][1] / N;
		components[k].reverse.w 	= EX[k][2] / N;
		W+=(components[k].bidir.w + components[k].forward.w + components[k].reverse.w  );
	}
	for (int k = 0; k < K; k++){
		components[k].bidir.w 		= components[k].bidir.w / W;
		components[k].forward.w 	= components[k].forward.w / W;
		components[k].reverse.w 	= components[k].reverse.w / W;
	}
}



double move_uniforom_support2(component * components, int K, int add, 
	segment * data, double move, double base_ll, int direction, double var, 
	vector<double> left, vector<double> right){
	//====================================
	random_device rd;
	mt19937 mt(rd());
	normal_distribution<double> step(0,var);
	//====================================
	int N 	 	= data->XN;
	double ll 	= nINF;
	double oldb, olda;
	for (int c = 0; c < K; c++){
		oldb  	= components[c].forward.b;
		olda 	= components[c].reverse.a;
		components[c].forward.b+=step(mt);
		components[c].forward.b=min(right[c],components[c].forward.b );
		if (components[c].forward.b < (components[c].bidir.mu + (1./components[c].bidir.l)  )  ){
			components[c].forward.b 	= components[c].bidir.mu + (2./components[c].bidir.l)  ;
		}
		ll 		= calc_log_likelihood(components, K+1, data  );
		if (ll < base_ll){
			components[c].forward.b 	= oldb;
		}else{
			base_ll 	= ll;
		}
		components[c].reverse.a+=step(mt);
		components[c].reverse.a = max(components[c].reverse.a, left[c]);
		ll 		= calc_log_likelihood(components, K+1, data  );
		if (ll < base_ll){
			components[c].reverse.a 	= olda;
		}else{
			base_ll 	= ll;
		}
		
	}	
	return base_ll;
			
}





double move_foot_print(segment * data, component * components, double base, int K, int add, double max_foot){
	random_device rd;
	mt19937 mt(rd());
	uniform_real_distribution<double> step(0,2);
	double current_ll;
	double threshold 	= pow(10,-5);
	double old_print, delta, base_ll, new_ll;
	double fvl, rvl, nfvl, nrvl;
	for (int c = 0; c < K; c++){
		old_print 	= components[c].bidir.foot_print;
		components[c].bidir.foot_print=0;
		delta 		= step(mt); 
		base_ll 	= 0, new_ll = 0;
		double start 	= components[c].bidir.mu - (1.0 / components[c].bidir.l) - components[c].bidir.si;
		double stop 	= components[c].bidir.mu + (1.0 / components[c].bidir.l) + components[c].bidir.si;
		int j 	= 0;
		while (j < data->XN and data->X[0][j]<start){
			j++;
		}
		int k=j;
		while (k < data->XN and data->X[0][j] < stop){
			k++;
		}

		for (int i = j; i < k; i++){
			fvl 	= components[c].bidir.pdf(data->X[0][i] -old_print , 1);
			rvl 	= components[c].bidir.pdf(data->X[0][i] +old_print , -1);
			
			nfvl 	= components[c].bidir.pdf(data->X[0][i]-delta, 1);
			nrvl 	= components[c].bidir.pdf(data->X[0][i]+delta, -1);
			
			base_ll+=LOG(fvl )*data->X[1][i];
			base_ll+=LOG(rvl )*data->X[2][i];				
			new_ll+=LOG(nfvl )*data->X[1][i];
			new_ll+=LOG(nrvl )*data->X[2][i];				
			
			
		}
		if (new_ll < base_ll){

			components[c].bidir.foot_print 	= old_print;
		}else{
			components[c].bidir.foot_print 	= delta;	
		}

	
	}
	double NLL 	= calc_log_likelihood(components, K+add, data);
	//cout<< NLL << " "<<base<<" "<<(NLL > base)<<endl;

	return base;
}

void update_elongation_supports(component * components, int K){
	for (int k = 0; k < K;k++)		{
		if (components[k].forward_neighbor!=NULL){
			components[k].forward.b 	= components[k].forward_neighbor->bidir.mu;
		}
		if (components[k].reverse_neighbor!=NULL){
			components[k].reverse.a 	= components[k].reverse_neighbor->bidir.mu;
		}
			
	}

}


//=========================================================
//For Classifier class / wrapper around EM

classifier::classifier(int k, double ct, int mi, double nm,
	double R_MU, double alpha_0, double beta_0,
	double alpha_1, double beta_1, double alpha_2,double alpha_3, double fp){
	foot_print 				= fp;
	K 						= k ;
	seed 					= true;
	convergence_threshold 	= ct;
	max_iterations 			= mi;
	noise_max 				= nm;
	p 						= 0.8;
	last_diff 				= 0;
	r_mu 					= R_MU;

	//=============================
	//hyperparameters
	ALPHA_0=alpha_0, BETA_0=beta_0, ALPHA_1=alpha_1, BETA_1=beta_1;
	ALPHA_2=alpha_2, ALPHA_3=alpha_3;

	move_l = true;

}
classifier::classifier(int k, double ct, int mi, double nm,
	double R_MU, double alpha_0, double beta_0,
	double alpha_1, double beta_1, double alpha_2,double alpha_3, bool MOVE, double fp){
	foot_print 				= fp;
	K 						= k ;
	seed 					= true;
	convergence_threshold 	= ct;
	max_iterations 			= mi;
	noise_max 				= nm;
	p 						= 0.8;
	last_diff 				= 0;
	r_mu 					= R_MU;

	//=============================
	//hyperparameters
	ALPHA_0=alpha_0, BETA_0=beta_0, ALPHA_1=alpha_1, BETA_1=beta_1;
	ALPHA_2=alpha_2, ALPHA_3=alpha_3;
	move_l 	= MOVE;
}
classifier::classifier(double ct, int mi, double nm,
	double R_MU, double alpha_0, double beta_0,
	double alpha_1, double beta_1, double alpha_2,double alpha_3, vector<vector<double>> IP, double fp){
	
	foot_print 				= fp;
	seed 					= true;
	convergence_threshold 	= ct;
	max_iterations 			= mi;
	noise_max 				= nm;
	p 						= 0.8;
	last_diff 				= 0;
	r_mu 					= R_MU;

	//=============================
	//hyperparameters
	ALPHA_0=alpha_0, BETA_0=beta_0, ALPHA_1=alpha_1, BETA_1=beta_1;
	ALPHA_2=alpha_2, ALPHA_3=alpha_3;
	init_parameters 		= IP;
	
}






classifier::classifier(){
	foot_print 	= 0;
};//empty constructor

int classifier::fit(segment * data, vector<double> mu_seeds ){
	//=========================================================================
	//for resets
	double fpp 		= foot_print;
	random_device rd;
	mt19937 mt(rd());
	uniform_real_distribution<double> dist_uni(data->minX,data->maxX);
	
	double l 	= data->maxX - data->minX;
	pi 	= sum(data->X[1], data->XN)/ data->N;
	double vl 	= 1.0 / l;
	double center 	= (data->start + data->stop) / 2.;
	
	if (K==0){
		//calc_likeihood coming from uniform model, only
		ll 	= 0;
		for (int i = 0; i < data->XN; i ++){
			if (pi > 0){
				ll+=(LOG(vl*(pi) )*data->X[1][i]);
				
			}
			if (pi < 1){
				ll+=(LOG(vl*(1-pi))*data->X[2][i]);
			}
		}
		converged=true;
		last_diff=0;
		components 	= new component[1];
		components[K].initialize(0., data, 0., 0. , noise_max, pi,0.);
		return 1;
	}
	int add 	= noise_max>0;
       
	components 	= new component[K+add];

	//initialize components
	for (int k = 0; k < K; k++){
		components[k].set_priors(ALPHA_0, BETA_0, ALPHA_1, BETA_1, ALPHA_2, ALPHA_3,data->N,K);
	}

	//===========================================================================
	//random seeds, initialize
	int i 	= 0;
	double mu;
	
	for (int k = 0; k < K; k++){
		if (mu_seeds.size()>0){
			i 	= sample_centers(mu_seeds ,  p);
			mu 	= mu_seeds[i];
			uniform_real_distribution<double> dist_uni2(mu-3,mu+3);
	
			mu 	= dist_uni2(mt);

		}else{
			mu 			= dist_uni(mt);
		}
		components[k].initialize(mu, data, K, data->SCALE , 0., 0.,foot_print);
		if (mu_seeds.size() > 0){
			mu_seeds.erase (mu_seeds.begin()+i);	
		}
	}
	
       
	if (add){
		components[K].initialize(0., data, 0., 0. , noise_max, pi, foot_print);
	}
 		
	//===========================================================================
	int t 			= 0; //EM loop ticker
	double prevll 	= nINF; //previous iterations log likelihood
	converged 		= false; //has the EM converged?
	double norm_forward, norm_reverse,N; //helper variables
	int u 			= 0;
	while (t < max_iterations && not converged){
		if (u > 100){
			check_mu_positions( components, K);
			u=0;
		}
		u++;
		ll 			= 0;
		//******
		//reset old sufficient statistics
		for (int k=0; k < K+add; k++){
			components[k].reset();
			if (components[k].EXIT){
				converged=false, ll=nINF;
				return 0;
			}
		       
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

			if (norm_forward > 0){
				ll+=LOG(norm_forward)*data->X[1][i];
			}
			if (norm_reverse > 0){
				ll+=LOG(norm_reverse)*data->X[2][i];
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
		if (abs(ll-prevll)<convergence_threshold){
			for (int k = 0; k < K; k++){
				double std 	= (components[k].bidir.si + (1.0 / components[k].bidir.l));
			}	
			converged=true;
		}
		if (not isfinite(ll)){
			ll 	= nINF;
			return 0;	
		}
		last_diff=abs(ll-prevll);

		prevll=ll;
		
		t++;
	}
	for (int k = 0; k < K; k++){
		double std 	= (components[k].bidir.si + (1.0 / components[k].bidir.l)) ;
		if (std > 10){
			ll 		= nINF;	
		}
	}	


	return 1;
}
vector<vector<double>> sort_mus(vector<vector<double>> X){
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

int classifier::fit_uniform_only(segment * data){
	//so in this case, there are a set of bidirectionals that are set and we are going to just try and maximize their 
	//elongation support first (and then maybe restimate parameters?)
	double foot_print 	=1;
	int K 			= init_parameters.size();
	init_parameters  = sort_mus(init_parameters);
	components 	= new component[K+1];
	//printf("--BEFORE--\n");
	vector<double> left;
	vector<double> right;
	double MU, L;
	left.push_back(data->minX);
	for (int k = 0; k < K; k++){
		MU 	= init_parameters[k][0], L 	= init_parameters[k][2];
		if (k > 0){
			right.push_back(MU);
		}
		if ( k < K-1 ){
			left.push_back(MU);
		}
	}
	right.push_back(data->maxX);

	for (int k =0; k < K;k++){
		components[k].initialize_with_parameters(init_parameters[k], data, K, left[k], right[k],foot_print);
		//components[k].print();
	}
			
	vector<double> empty;
	components[K].initialize_with_parameters(empty, data, K, data->minX,data->maxX,foot_print);

	bool converged=false;
	int t=0;
	ll 	= calc_log_likelihood(components, K+1, data  );
	double prevll 	= nINF;
	int times 		= 0;
	int changes 	= 0;
	double var 		= 10;

	while (not converged and t < max_iterations){
		
		update_weights_only(components, data, K, 1);
		ll 	= calc_log_likelihood(components, K+1, data  );
		ll 	= move_uniforom_support2( components, K, 1, data, move,  ll, 1, var, left, right);
		//update_weights_only(components, data, K, 1);
		//ll 	= calc_log_likelihood(components, K+1, data  );
		
		//printf("%f\n",ll );
		if (abs(ll-prevll) < convergence_threshold){
			changes++;
		}else{
			changes=0;
		}
		if (changes >20){
			converged=true;
		}

		prevll=ll;

		t++;
	}
	//now we can probably remaximize...but maybe save that for a later date...


	return 1;

}
int classifier::fit_uniform_only2(segment * data){
	
	
	return 1;

}
string classifier::print_out_components(){
	string text 	= "";
	for (int k = 0; k < K; k++){
		text+=components[k].write_out();
	}
	return text;
}
void print_components(component * components, int K, double ll){
	printf("-------------------- === %f\n", ll);
	for (int c = 0; c < K; c++){
		components[c].print();
	}
}

void sort_vector(double X[], int N){
	bool sorted=true;
	while (sorted){
		sorted=false;
		for (int i = 0; i < N-1; i++){
			if (X[i] > X[i+1]){
				double cp = X[i];
				X[i] 		= X[i+1];
				X[i+1] 		= cp;
				sorted 		= true;
			}
		}
	}
}

int classifier::fit2(segment * data, vector<double> mu_seeds, int topology,
	 int elon_move ){


	if (K == 0){
		ll 			= 0;
		double 	l 	= (data->maxX-data->minX);
		double pos 	= 0;
		double neg 	= 0;
		for (int i = 0; i < data->XN; i++){
			pos+=data->X[1][i];
			neg+=data->X[2][i];
		}
		double pi 	= pos / (pos + neg);
		for (int i = 0; i < data->XN; i++){
			if (pi > 0){
				ll+=log(pi / l)*data->X[1][i];
			}
			if (pi < 1){
				ll+=log((1-pi) / l)*data->X[2][i];		
			}
		}

		components 	= new component[1];

		return 1;
	}
	//=========================================================================
	//for resets
	random_device rd;
	mt19937 mt(rd());
	uniform_real_distribution<double> dist_uni(data->minX,data->maxX);
	uniform_real_distribution<double> dist_sample_temp(0,1);
	
	int add 	= noise_max>0;
	components 	= new component[K+add];

	//initialize components
	for (int k = 0; k < K; k++){
		components[k].set_priors(ALPHA_0, BETA_0, ALPHA_1, BETA_1, ALPHA_2, ALPHA_3,data->N, K);
	}

	//===========================================================================
	//random seeds, initialize
	int i 	= 0;
	double mu;
	//get mus to order
	double mus[K];
	for (int k = 0; k < K; k++){
		double U 	= dist_sample_temp(mt);
		if (mu_seeds.size()>0 and U < 1.0  ){
			i 	= sample_centers(mu_seeds ,  p);
			mu 	= mu_seeds[i];
			if (r_mu > 0){
				normal_distribution<double> dist_r_mu(mu, r_mu);
				mu 		= dist_r_mu(mt);
			}
		}else{
			mu 			= dist_uni(mt);
		}
		
		mus[k] 	= mu;
		if (mu_seeds.size() > 0 and U < 0.5){
			mu_seeds.erase (mu_seeds.begin()+i);	
		}
	}
	sort_vector(mus, K);
	//sort_vector(mus,K);
	double reverse_bound=data->minX, forward_bound=data->maxX;
	int pad 	= 500;
	for (int k = 0; k < K;k++){
		if (k==0){
			reverse_bound=data->minX;
		}else if(k>0){
			reverse_bound=mus[k-1]+(pad/data->SCALE);
		}
		if (k+1 < K){
			forward_bound=mus[k+1]-(pad/data->SCALE);
		}else{
			forward_bound=data->maxX;
		}

		components[k].initialize_bounds(mus[k], data, K, data->SCALE , 0., topology,foot_print, forward_bound, reverse_bound);
		
	}
	sort_components(components, K);
	
	for (int k = 0; k < K; k++){
		if (k > 0){
			components[k].reverse_neighbor 	= &components[k-1];
		}else{
			components[k].reverse_neighbor 	= NULL;
		}
		if (k+1 < K){
			components[k].forward_neighbor 	= &components[k+1];
		}else{
			components[k].forward_neighbor 	= NULL;
		}
	}

       
	if (add){
		components[K].initialize_bounds(0., data, 0., 0. , noise_max, pi, foot_print, data->minX, data->maxX);
	}
	if (topology<2 and elon_move ){
		update_l(components, data, K);
	}
	//===========================================================================
	int t 			= 0; //EM loop ticker
	double prevll 	= nINF; //previous iterations log likelihood
	converged 		= false; //has the EM converged?
	double norm_forward, norm_reverse,N; //helper variables
	vector<double> left;
	vector<double> right;
	double MU, L;
	int u 			= 0;
	while (t < max_iterations && not converged){
		//******
		//reset old sufficient statistics

		for (int k=0; k < K+add; k++){
			components[k].reset();
			if (components[k].EXIT){
				converged=false, ll=nINF;
				return 0;
			}
		       
		}
		
		//******
		//E-step
		ll 	= 0;
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
			if (norm_forward > 0){
				ll+=LOG(norm_forward)*data->X[1][i];
			}
			if (norm_reverse > 0){
				ll+=LOG(norm_reverse)*data->X[2][i];
			}
			//now we need to add the sufficient statistics
			for (int k=0; k < K+add; k++){
				if (norm_forward){
					components[k].add_stats(data->X[0][i], data->X[1][i], 1, norm_forward);
				}
				if (norm_reverse){
					components[k].add_stats(data->X[0][i], data->X[2][i], -1, norm_reverse);
				}
				//components[k].print();
			}
		}

		//******
		//M-step
		//get normalizing constant
		N=0;
		for (int k = 0; k < K+add; k++){
			N+=(components[k].get_all_repo());
		}
		
		for (int k = 0; k < K+add; k++){
			components[k].update_parameters(N, K);
		}
		
		if (abs(ll-prevll)<convergence_threshold){
			converged=true;
		}
		if (not isfinite(ll)){
			ll 	= nINF;
			return 0;	
		}
		if (u > 100 ){
			sort_components(components, K);
			check_mu_positions(components, K);
			if (elon_move){
				update_j_k(components,data, K, N);
				update_l(components,  data, K);
			}
			u 	= 0;
		}

		u++;
		
		last_diff=abs(ll-prevll);

		prevll=ll;
		
		t++;
	}

	return 1;
}














