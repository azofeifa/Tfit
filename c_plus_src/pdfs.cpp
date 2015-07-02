#include <math.h>
#include <algorithm>
#include <iostream>
#include "pdfs.h"
using namespace std;
double IN(double x){
	return exp(-pow(x,2)/2.)/(sqrt(2*M_PI));
}
double IC(double x){
	return 0.5*(1+erf(x/sqrt(2)));
}
double R(double x){
	double N = IC(x);
	double D = IN(x);
	if (D < pow(10,-25)){
		return 1.0 / pow(10,-25);
	}
	return (1.0 - N) / D;
}


//====================================================================
//Exponentially Modiefied Gaussian
EGU::EGU(){};
void EGU::set(double MU, double S, double L, double B, double W, double PI_EMG, double PI_UNI ){
	//weights
	mu 		= MU;
	s 		= S;
	l  		= L;
	b 		= B;
	w 		= W;
	pi_emg 	= PI_EMG;
	pi_uni 	= PI_UNI;

}
double EGU::pdf_U(double x, int f){
	if (mu <= x and x<= b){
		double vl 	= (1.- w ) / (b-mu);
		if (f 	== 0){
			return vl * pi_uni;
		}
		return vl * (1-pi_uni);
	}
	return 0.;
		
}


double EGU::pdf_EG(double x, int f){
	if (f==0){
		return w*pi_emg*(l/2.)*exp((l/2.)*(2*mu + l*pow(s,2) - 2*x))*erfc((mu+ l*pow(s,2)-x)/(sqrt(2)*s) );
	}else{
		return w*(1-pi_emg)*(l/2.)*exp((l/2.)*(2*x - 2*mu + l*pow(s,2) ))*erfc((x - mu+ l*pow(s,2))/(sqrt(2)*s) );	
	}
}

double EGU::cond_y(double z, int f){
	if (f==0){
		return max(0., z - mu - l*pow(s,2) +  (s/ R(l*s - ((z-mu)/s) ))) ;
	}
	return max(0. , (s/ R(l*s + ((z-mu)/s) ) ) - (l*pow(s,2) -mu+z)) ;
}
double EGU::cond_y_sq(double z, int f){
	if (f==0){
		return max(0. , pow(l,2)*pow(s,4) + pow(s,2)*(2*l*(mu-z) + 1) + pow(mu-z,2) - ((s*(l*pow(s,2) + mu - z))   / R(l*s - ((z-mu)/s)) ));
	}
	return max(0., pow(l,2)*pow(s,4) + pow(s,2)*(2*l*(z-mu) + 1) + pow(mu-z,2) + ((s*(l*pow(s,2) - mu + z))   / R(l*s + ((z-mu)/s) )) );
}

double EGU::cond_x(double z, int f){
	if (f==0){
		return z - cond_y(z, f);
	}else{
		return z + cond_y(z, f);
	}
}

double EGU::cond_x_sq(double z, int f ){
	return pow(cond_x(z,f), 2) + cond_y_sq(z,f) - pow(cond_y(z,f),2);
}
			

void EGU::set_partial(double MU, double S, double L, double PI_EMG, double PI_UNI){
	mu 		= MU;
	s 		= S;
	l 		= L;
	pi_emg 	= PI_EMG;
	pi_uni 	= PI_UNI;
}
void EGU::set_weights(double W){
	//weights
	w 	= W;
}

//====================================================================
//====================================================================
//====================================================================
//Double Geometric
DGU::DGU(){};
void DGU::set( double I, double U, double S, double L, double W){
	i=I;
	u=U;
	s=S;
	l=L;
	w=W;
}

double DGU::pdf_U(double x){
	if (l >0 and 0<=x and x<=l){
		return ((1.- w)/ l);
	}else if(l <0 and l<=x and x<=0){
		return ((1.- w)/abs(l));
	}
	return 0.;
}
double DGU::pdf_DG(double k){	
	return w*u*pow(1-u, max(0., -k))*s*pow(1-s, max(0.,k))/ (u+s-(u*s));
}








