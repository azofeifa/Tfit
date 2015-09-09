#include <math.h> 
#include <limits>
#include "model_single.h"
#include "load.h"
#include "template_matching.h"
#include <iostream>
#include <algorithm>
#include <unistd.h>
#include <random>

using namespace std;
NORM::NORM(){}
NORM::NORM(double MU, double SI, double W){
	mu=MU, si=SI, w=W;
}
string NORM::print(){
	return "N, mu:" + to_string(mu) + ",si:" + to_string(si) +",w:"+to_string(w);
}
double NORM::pdf(double x){
	return (w/(si*sqrt(2*M_PI)))*exp(-pow(x-mu,2)/(2*pow(si,2) ));
}

ELON::ELON(){}
ELON::ELON(double A, double B, double W){
	a=A, b=B, w=W;
}
double ELON::pdf(double x){
	if (w == 0){
		return 0;
	}
	if (x >= a and x<= b){
		return (w / (b-a));
	}
}
string ELON::print(){
	return "U, a:" + to_string(a) + ",si:" + to_string(b) +",w:"+to_string(w);

}

void NLR::print(){
	string out = "";
	out+=loading.print()+ "\n";
	out+=forward.print()+ "\n";
	out+=reverse.print()+ "\n";
	cout<<out;

}

NLR::NLR(){
	 EX=0, EX2=0, WN=0, WR=0, WL=0;
}
void NLR::init(int type, int K, segment * data, double scale){
	if (type==0){
		wn 	=0, wl=0, wr=1.0;
	}else if(type==1){
		wn 	=1./double(2*K), wl=0, wr=1./double(2*K);
	}else if(type==2){
		wn 	=1./double(2*K), wl=1./double(2*K), wr=0;
	}else{
		wn 	=1./double(3*K), wl=1./double(3*K), wr=1./double(3*K);	
	}
	random_device rd;
	mt19937 mt(rd());
	uniform_real_distribution<double> dist_mu(data->minX, data->maxX);
	uniform_real_distribution<double> dist_si(500, 5000);
	
	mu 	= (data->maxX + data->minX) / 2.;
	si 	= 10;
	l 	= data->minX;
	r 	= data->maxX;
	loading 	= NORM(mu, si, wn);
	forward 	= ELON(mu,r, wr);
	reverse 	= ELON(l,mu, wl);


}
double NLR::pdf(double x){
	return loading.pdf(x) + forward.pdf(x) + reverse.pdf(x);
}

void NLR::addSS(double x , double y, double norm){
	double lp 	= loading.pdf(x)/norm;
	printf("x: %f, EX: %f, p: %f\n", x, lp, loading.pdf(x) );
	double fp 	= forward.pdf(x)/norm;
	double rp 	= reverse.pdf(x)/norm;
	EX+=x*y*lp;
	EX2+=pow(x-mu,2)*y*lp;
	WN+=y*lp;

	WR+=y*fp;
	WL+=y*rp;
}
void NLR::resetSS(){
	 EX=0, EX2=0, WN=0, WR=0, WL=0;	
}
double NLR::get_all(){
	return WN+WR+WL;
}
void NLR::set_new_parameters(double N){
	wn 			= WN / N, wl 	= WL / N, wr 	= WR / N;
	mu  		= EX / WN;
	si 			= sqrt((EX2+0.01) /(WN+0.01));
	loading 	= NORM(mu, si, wn);
	forward 	= ELON(mu,r, wr);
	reverse 	= ELON(l,mu, wl);

}

classifier_single::classifier_single(double ct, int mi, int k, int TY, double SC){
	covergence_threshold=ct, max_iterations=mi;
	K = k, type = TY;
	scale 	= SC;
}

classifier_single::classifier_single(){}


double calc_loglikelihood(segment * data, NLR * components, int K){
	double ll 	= 0;
	double S;
	for (int i =0; i < data->XN; i++){
		S 	= 0;
		for (int k = 0; k < K; k++){
			S+=(components[k].pdf(data->X[0][i]));
		}
		ll+=log(S)*data->X[1][i];
	}
	return ll;
}

double classifier_single::fit(segment * data){
	components 	= new NLR[K];
	int t 			= 0;
	double prevll 	= nINF;
	if (type == 0){
		double vl 			= 1.0 / (data->maxX - data->minX);
		ll 			= 0;
		for (int i = 0; i < data->XN;i++){
			ll+=log(vl)*data->X[1][i];
		}
		return ll;
	}
	//random initialization of parameters
	for (int k = 0; k < K; k++ ){
		components[k].init(type, K, data, scale);
	//	components[k].print();
	}
	ll 				= calc_loglikelihood(data, components, K);
	bool converged 	= false;
	double norm, N;

	printf("\n---------------\n");
	while (t < max_iterations and not converged){
		if (not isfinite(ll)){
			ll 		= nINF;
			return ll;
		}
		for (int k = 0; k < K; k++ ){
			components[k].resetSS();
			components[k].print();
		}	
		for (int i = 0 ; i < data->XN; i++){
			norm 	= 0;
			for (int k = 0; k < K;k++){
				norm+=components[k].pdf(data->X[0][i] );
			}
			if (norm > 0){
				for (int k = 0; k < K;k++){
					components[k].addSS(data->X[0][i], data->X[1][i], norm);
				}
			}

		}
		N 	= 0;
		for (int k = 0;k < K;k++){
			N+=components[k].get_all();
		}
		for (int k = 0;k < K;k++){
			components[k].set_new_parameters(N);
		}
		ll 				= calc_loglikelihood(data, components, K);
		if (abs(ll-prevll) < covergence_threshold){
			for (int c = 0; c < K; c++){
				components[c].print();
				if (components[c].si > (2000/scale)){
					ll 	=nINF;
				}
			}
			converged=true;
		}
		prevll 			= ll;
		t++;
	}

	return ll;
}


