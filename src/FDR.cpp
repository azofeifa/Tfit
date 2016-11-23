#include "FDR.h"
#include "load.h"
#include "model.h"
#include <limits>
#include <iostream>
#include <algorithm>
#include <random>
#include "omp.h"
#include <cmath>
#include <vector>
#include <string>
#include "BIC.h"
#include <math.h> 
using namespace std;

normal::normal(){}
normal::normal(double X, double MU, double STD){
  mean=MU, x=X, std=STD, c2 = sqrt(2);
}
double normal::cdf(double x){
  double pv = 0.5;
  double Z  =erf( (x-this->mean) / (this->std*c2));
  pv        = 0.5*(1+Z);
  return pv;
}
double normal::pdf(double x){
  return 1.0 / (sqrt(2*M_PI)*this->std  )*exp(-pow(x-this->mean,2)/(2*pow(this->std,2) ));  
}

exponential::exponential(){};
exponential::exponential(double MU, double LAMBDA){
  this->mu = MU , this->lambda = LAMBDA;
}
double exponential::pdf(double x){
  if (x>=mu){
    return this->lambda*exp(-this->lambda*(x-mu) ); 
  }
  return 0.0;
}

pareto::pareto(){};
pareto::pareto(double MU, double ALPHA){
  this->mu = MU , this->alpha=ALPHA;
};

double pareto::pdf(double x){
  if (x>=mu){
    return this->alpha*pow(this->mu,alpha) / pow(x, this->alpha+1);
  }
  return 0.0;
}

bool check_value(float x){
  if (isnan(x) ){
    return false;
  }
  if (not isfinite(x) ){
    return false;
  }
  return true;
}

bool EM(vector<vector<double>> X , double & mean, double & sigma, double & w, double & c , bool EXP ){
  mean = 0.3 , sigma = 0.1, w = 0.5;
  if (EXP){
    c=0.5;
  }else{
    c=1.34;
  }
  normal N(1.0, mean,sigma) ;
  exponential E(mean, c);
  pareto P(mean,c);
  
  int t=0, T=300;
  double prevw = -0.1;
  bool converged=false;
  while (t < T and not converged  ){
    double EX=0.0,EX2=0.0, EY=0.0,R1=0.0, R2=0.0,p1=0.0,p2=0.0,r1=0.0,r2=0.0;
    for (int i = 0 ; i < X.size();i++){
      
      double x = X[i][0], y = X[i][1];
      p1= (1-w)*N.pdf(x);
      if (EXP){
	p2=w*E.pdf(x);
      }else{
	p2=w*P.pdf(x);
      }
      if ((p1 + p2) > 0){
	r1= p1 / (p1 + p2) ,r2=p2 / (p1 + p2);
	R1 += r1*y,R2+=r2*y;
	EX +=(r1*x)*y;
	if (EXP){
	  EY +=(r2*x)*y;
	}else{
	  EY += (log(x) - log(N.mean))*r2*y;
	}
	EX2 +=(r1*pow(x-N.mean,2))*y;
      }
    }
    
    N.mean = EX /R1;
    N.std = sqrt(EX2 / R1);
    if (EXP){
      E.lambda = R2/EY ;
    }else{
      P.alpha = R2/EY;
    }
    if( not check_value(N.mean)  or not check_value(N.std)  ){
      converged = false;
      break;
    }
    w = R2 / (R1 + R2);
    if (abs(w-prevw) < pow(10,-6)){
      converged=true;
    } 
    prevw=w, t+=1;
  }
  if(not converged){
    mean = 0.3 , sigma    = 0.05;
  }else{
    
    mean = N.mean , sigma = N.std*15;
  }
  return converged;
} 



slice_ratio::slice_ratio(){};

slice_ratio::slice_ratio(double ST, double SP, int BINS){
  this->start = ST, this->stop = SP, this->bins=BINS;
  double step = (stop - start)/double(this->bins);
  for (double c = this->start; c < this->stop ; c+=step){
    vector<double> row  = {c, 0.0}; //step, N, X, X^2
    this->XY.push_back(row);
  }
};
int slice_ratio::get_closest(double x){
  int c = 0;
  while (c < this->XY.size() and this->XY[c][0] < x ){
    c+=1;
  }
  if (c > 0){
    return c -1;
  }
  return 0;

}
void slice_ratio::insert(double y){
  int c = this->get_closest(y);
  this->XY[c][1]+=1;
}




void slice_ratio::set(double pval){
  this->mean=0.1, this->std=0.1, this->w=0.5, this->c=1.0;
  this->converged = EM(this->XY, this->mean, this->std, this->w, this->c, 1);
  this->norm_all  = normal(0.0, this->mean, this->std);
  double z        = this->norm_all.mean - this->norm_all.std;
  double pv       = this->norm_all.cdf(z);
  while (pv < (1.0 - pval) and z < 100.0){
    z+=0.01;
    pv=this->norm_all.cdf(z);
  }
  threshold=z;
}

void slice_ratio::set_2(double pval){
  this->norm_all  = normal(0.0, this->mean, this->std);
  double z        = this->norm_all.mean - this->norm_all.std;
  double pv       = this->norm_all.cdf(z);
  while (pv < (1.0 - pval) and z < 100.0){
    z+=0.01;
    pv=this->norm_all.cdf(z);
  }
  threshold=z;
}

double slice_ratio::pvalue(double y){
  double pv = 1.0-this->norm_all.cdf(y);
  return pv;
}

slice_ratio get_slice(vector<segment *> segments, int N, double CC, params * P){

  double sigma, lambda, fp, pi, w, window, pval_threshold,ns;
 
  window        = stod(P->p["-pad"]), ns=stod(P->p["-ns"]) ;
  sigma         = stod(P->p["-sigma"])/ns , lambda= ns/stod(P->p["-lambda"]);
  fp            = stod(P->p["-foot_print"])/ns , pi= stod(P->p["-pi"]), w= stod(P->p["-w"]);
  pval_threshold= stod(P->p["-bct"]) ;
  int CN     = segments.size();
  double min_x  = -1 , max_x = -1, n = 0;
  random_device rd;
  mt19937 mt(rd());
  default_random_engine generator;
  uniform_real_distribution<double> distribution(0,1);
  vector<double> XY(N);
  vector<double> CovN(N);
  for (int i = 0 ; i < XY.size(); i++){
    XY[i]=0.0, CovN[i]=0.0;
  }
  #pragma omp parallel for
  for (int n = 0 ; n < N ; n++){
    double U       = distribution(mt);
    double U2      = distribution(mt);
    int NN         = int(U*(CN-1));
    segment * data = segments[NN];
    int c          = U2*int(data->XN);
    int j = c,  k  = c;
    double N_pos = 0 , N_neg =0 ;
    while (j > 0 and (data->X[0][c] - data->X[0][j] )< window){
      N_pos+=data->X[1][j];
      N_neg+=data->X[2][j];
      j--;
    }
    while (k < data->XN and (data->X[0][k] - data->X[0][c] )< window  ){
      N_pos+=data->X[1][k];
      N_neg+=data->X[2][k];
      k++;
    }
    CovN[n] = N_pos + N_neg;
    if (N_pos + N_neg > CC and (data->X[0][k] - data->X[0][j]) > 1.75*window  ){
      
      double val =  BIC3(data->X,  j,  k,  c, N_pos,  N_neg, sigma , lambda, fp , pi, w);
      if (val >0 ){
        XY[n]=val;
      }
    }
  }
  double S1=0.0, S2=0.0, TT=0.0;
  for (int n = 0 ; n < CovN.size();n++){
    if (CovN[n]>0){
      S1+=CovN[n], S2+=pow(CovN[n],2), TT+=1; 
    }
  }
  //---------------
  for (int n = 0 ; n < XY.size(); n++){
    if (min_x < 0 or XY[n] < min_x ){
      min_x = XY[n];
    }
    if (max_x < 0 or XY[n] > max_x ){
      max_x = XY[n];
    }
  }
  //-------------
  slice_ratio SC(min_x,max_x,400);
  for (int n = 0 ; n < XY.size();n++){

    if (XY[n] > 0.0){
      SC.insert(XY[n] );
    }
  }
  SC.set(pval_threshold);
  return SC;
}
