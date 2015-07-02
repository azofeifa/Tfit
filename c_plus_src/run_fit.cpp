#include "load.h"
#include "math.h"
#include "pdfs.h"
#include <algorithm>
#include <iostream>
#include <chrono>
using namespace std;
double getMean(double ** X, int N, int s){
	double S 	= 0;
	double n 	= 0;
	for (int i = 0; i < N; i ++){
		S+=X[i][0]*X[i][s];
		n+=X[i][s];
	}
	return S/n;
}
double getStd(double ** X,  int N, int s){
	double S 	= 0;
	double n 	= 0;
	double mean = getMean(X,N,s);
	for (int i = 0; i < N; i ++){
		S+=pow(X[i][0]-mean,2)*X[i][s];
		n+=X[i][s];
	}
	return sqrt(S/n);	
}
double getWeightedMean(double * X,double * Y, int N){
	double S 	= 0;
	double n 	= 0;
	for (int i = 0; i < N; i ++){
		if (not std::isnan(X[i]) and not std::isnan(Y[i]) and not std::isinf(X[i]) and not std::isinf(Y[i]) ){
			S+=(X[i]*Y[i]);
			n+=Y[i];
		}

		
	}
	return S/n;
}
double getWeightedSum(double *X, double *Y, int N){
	double S 	= 0;
	for (int i = 0; i < N; i++){

		if (not std::isnan(X[i]) and not std::isnan(Y[i])){
			S+=X[i]*Y[i];
		}
	}
	return S;
}
double getSum(double * X, int N){
	double S 	= 0;
	for (int i = 0; i < N; i++){
		if (not std::isnan(X[i])){
			S+=X[i];
		}
	}
	return S;
}

double compute_E_variance(double * X, double * X_sq, double * Y, double mu, int N){
	double s 	= 0;
	double sq 	= 0;
	double n 	= 0;
	for (int i = 0; i < N; i ++){
		if (not std::isnan(X[i]) and not std::isnan(X_sq[i]) and not std::isnan(Y[i]) ){
			s+=X[i]*Y[i];
			sq+=X_sq[i]*Y[i];
			n+=Y[i];
		}
	}
	return sqrt((sq - 2*mu*s + n*pow(mu,2))/ n);
}




bool fit_DG(annotation * data){
	int N 	= data->N;
	
	double r[N][2];  //responsibility terms, 0 is DG, 1 is U
	double minX 	= data->D[0][0];
	double maxX 	= data->D[N-1][0];
	int s;
	if (data->strand == "+"){
		s 	= 0;
	}else{
		s 	= 1;
	}

	int st,sp, i;
	double B;

	

	if (s == 0){
		i 	= 1;
		B 	= maxX;
		sp 	= N/4;
		st 	= 0;
	}else{
		i 	= sp-1;
		B 	= minX;
		sp 	= N;
		st 	= N*0.75;
	}
	double ** D = data->D;
	double I, DG_val, U_val,norm,alpha,beta, U,W,S ;
	double NN;
	double qj_norm ;
	double L,ll, prevll;
	double bestll 	= 0;

	int t;
	bool converged 	= false;
	DGU * bestrv 	= new DGU;
	while (i < sp and i > st){
		I 	= D[i][0];
		U 		= 0.001;
		W 		= 0.5;
		S 		= 0.001;

		if (s==0){
			L 	= maxX - I;
		}else{
			L 	= minX - I;
		}
		//=================================
		//compute responsibilities
		DGU rv;

		rv.set(I, U, S, L, W); 	
		prevll = 0;
		t=0;
		converged 	= false;
		while (not converged and t< 500){
			t++;
			qj_norm = 0;
			NN 		= 0;
			alpha 	= 0;
			beta 	= 0;
			ll 		= 0;
			for (int j = 0; j < N; j ++){
				DG_val 	= rv.pdf_DG(D[j][0]-I);
				U_val 	= rv.pdf_U(D[j][0]-I);
				ll 		+=log(DG_val + U_val)*D[j][s+1];
				norm 	= DG_val + U_val;
				if (norm){
					r[j][0] 	= DG_val/norm;
					r[j][1] 	= U_val/norm;
				}else{
					r[j][0] 	= 0;
					r[j][1] 	= 0;
				}
				qj_norm+=(r[j][0]*D[j][s+1]);
				NN+=D[j][s+1];
			}

			//=================================
			//compute alpha and beta terms
			for (int j =0; j < N;j++){
				if (D[j][0]<=I){
					alpha+=(I-D[j][0])*(r[j][0]/qj_norm)*D[j][s+1];
				}
				else if (D[j][0]>I){
					beta+=(D[j][0]-I)*(r[j][0]/qj_norm)*D[j][s+1];
				}
			}
			U 			= 2./(1+2*alpha + sqrt(1+4*alpha*beta));
			S 			= 2./(1+2*beta + sqrt(1+4*alpha*beta));
			W 			= qj_norm/(NN);
			
			rv.set(I, U, S, L, W); 	
			if (abs((prevll-ll)) < 0.0001){
				converged 	= true;
			}
			prevll 		= ll;
		}
		if (bestll == 0 or ll > bestll){
			bestll 	= ll;
			bestrv->set(I,U,S,L,W);
		}
		if (s==0){
			i++;
		}else{
			i--;
		}
	}
	data->rv_dgu_best 	= bestrv;
	return 1;
}
bool fit_EG(annotation * data, int K){
	//==========================================
	//alg parameters
	int max_iter 	= 100;
	int t 			= 0;
	int i,k;
	//==========================================
	//misc parameters
	bool converged 	= false;
	double z,yf, yr, norm,tot;
	double minz 	= data->D[0][0];
	double maxz 	= data->D[data->N-1][0];
	string strand 	= data->strand;
	double sum_EMG_f, sum_EMG_r, sum_U_f,sum_U_r, PI_EMG, PI_UNI;
	double  Mu, Lambda, Sigma;
	//==========================================
	//instantiate some arrays for E-step calcs 
	EGU rvs[K];
	double mixture_E[2][K*2][data->N];
	double EZ[2][K][3][data->N];


	//==========================================
	//initialize random parameter values for EGU
	double var, init_mu, init_lam;
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	
	std::default_random_engine generator(seed);
  	for (k =0; k < K; k++){
		if (strand == "+"){
			var 		= 50.;
			init_lam 	= 0.001;
			std::uniform_real_distribution<double> n_dist(minz, maxz);
			init_mu 	= n_dist(generator);

			rvs[k].set(init_mu, var, init_lam, maxz, 0.5, 0.5,0.5 );
		}else{
			var 		= 50.;
			init_lam 	= 0.001;
			std::uniform_real_distribution<double> n_dist(minz, maxz);	
			init_mu 	= n_dist(generator);
			rvs[k].set(init_mu, var, init_lam, maxz, 0.5, 0.5,0.5 );
		}
	}


	//==========================================
	//start
	while (t< max_iter and not converged){
		//===============
		//E-step:
		for (i=0; i < data->N; i++){
			z 	= data->D[i][0];
			yf 	= data->D[i][1];
			yr 	= data->D[i][2];
			if (yf){//there is a forward strand data point so compute 
				norm= 0;
				for (k=0; k < K;k++){
					//this is the E-step for component weights
					mixture_E[0][k*2][i] 		= rvs[k].pdf_EG(z, 0); //EM
					mixture_E[0][k*2+1][i] 		= rvs[k].pdf_U(z, 0); //forward uniform
					norm+=(mixture_E[0][k*2][i] + mixture_E[0][k*2+1][i]);
					//this is the E-step for the EMG
					if (mixture_E[0][k*2][i]){//only compute if we need to
						EZ[0][k][0][i] 			= rvs[k].cond_x(z, 0);
						EZ[0][k][1][i]			= rvs[k].cond_x_sq(z, 0);
						EZ[0][k][2][i] 			= rvs[k].cond_y(z, 0);
					}else{//sense we didn't compute this, we need to flag it, so we don't sum it up later...
						EZ[0][k][0][i] 			= numeric_limits<double>::quiet_NaN();
						EZ[0][k][1][i] 			= numeric_limits<double>::quiet_NaN();
						EZ[0][k][2][i] 			= numeric_limits<double>::quiet_NaN();							
					}
				}
				//renormalize for numerical stability
				for (k=0; k < K; k++){
					if (norm){
						mixture_E[0][k*2][i] 	/=norm;
						mixture_E[0][k*2+1][i]	/=norm;
						mixture_E[0][k*2][i] 	*=yf;
						mixture_E[0][k*2+1][i] 	*=yf;

					}else{
						mixture_E[0][k*2][i]	=0;
						mixture_E[0][k*2+1][i] 	=0;
					}
				}
			}else{
				for (k=0; k < K;k++){
					mixture_E[0][k*2][i] 	= 0;
					mixture_E[0][k*2+1][i]	= 0;
					EZ[0][k][0][i] 			= 0;
					EZ[0][k][1][i] 			= 0;
					EZ[0][k][2][i] 			= 0;
				}
			}
			if (yr){//there is a reverse strand data point so compute
				norm 	= 0;
				for (k=0; k < K;k++){
					//this is the E-step for component weights
					mixture_E[1][k*2][i] 		= rvs[k].pdf_EG(z, 1); //EM 
					mixture_E[1][k*2+1][i] 		= rvs[k].pdf_U(z, 1); //forward uniform
					norm+=(mixture_E[1][k*2][i] + mixture_E[1][k*2+1][i]);
					//this is the E-step for the EMG
					if (mixture_E[1][k*2][i]){
						EZ[1][k][0][i] 			= rvs[k].cond_x(z, 1);
						EZ[1][k][1][i] 			= rvs[k].cond_x_sq(z, 1);
						EZ[1][k][2][i] 			= rvs[k].cond_y(z, 1);
					}else{//sense we didn't compute this, we need to flag it, so we don't sum it up later...
						EZ[1][k][0][i] 			= numeric_limits<double>::quiet_NaN();
						EZ[1][k][1][i] 			= numeric_limits<double>::quiet_NaN();
						EZ[1][k][2][i] 			= numeric_limits<double>::quiet_NaN();			
					}
				}
				//renormalize for numerical stability
				for (k=0; k < K; k++){
					if (norm){
						mixture_E[1][k*2][i] 	/=norm;
						mixture_E[1][k*2+1][i]	/=norm;
						mixture_E[1][k*2][i] 	*=yr;
						mixture_E[1][k*2+1][i] 	*=yr;
					
					}else{
						mixture_E[1][k*2][i]	=0;
						mixture_E[1][k*2+1][i] 	=0;
					}
				}
			}else{
				for (k=0; k < K;k++){
					mixture_E[1][k*2][i] 	= 0;
					mixture_E[1][k*2+1][i]	= 0;
					EZ[1][k][0][i] 			= 0;
					EZ[1][k][1][i] 			= 0;
					EZ[1][k][2][i] 			= 0;
				}
			}
		}//phew...that is a big E-step
		//===================
		//M-step:
		tot=0;
		for (k=0;k<K;k++){
			//double W_p, W_ef, W_er, Pi_p, Pi_ef, Pi_er, Mu, Sigma, Lambda, L_f, L_r;
			sum_EMG_f 		= getSum(mixture_E[0][k*2], data->N);
			sum_EMG_r 		= getSum(mixture_E[1][k*2], data->N);
			sum_U_f 		= getSum(mixture_E[0][k*2+1], data->N);
			sum_U_r 		= getSum(mixture_E[1][k*2+1], data->N);
			tot+=(sum_EMG_f + sum_EMG_r + sum_U_f 	+ sum_U_r );
			PI_EMG 			= sum_EMG_f/(sum_EMG_f 	+ sum_EMG_r);
			PI_UNI 			= sum_U_f/(sum_U_f 		+ sum_U_r);
			//============================================================
			//Set EMG Parameters
			Mu 				= getWeightedMean(EZ[0][k][0],mixture_E[0][k*2], data->N)*PI_EMG + getWeightedMean(EZ[1][k][0],mixture_E[1][k*2], data->N)*(1. - PI_EMG);
			Lambda 			= min(abs(1.0 /(getWeightedMean(EZ[0][k][2],mixture_E[0][k*2], data->N)*PI_EMG + getWeightedMean(EZ[1][k][2],mixture_E[1][k*2], data->N)*(1. - PI_EMG))), 3.);
			Sigma 			= compute_E_variance(EZ[0][k][0], EZ[0][k][1], mixture_E[0][k*2], Mu, data->N);//*PI_EMG + compute_E_variance(EZ[1][k][0], EZ[1][k][1], mixture_E[1][k*2], Mu, data->N)*(1.-PI_EMG) ;
			cout<<PI_EMG<<endl;
			rvs[k].set_partial(Mu, Sigma, Lambda,PI_EMG, PI_UNI);
		}
		//============================================================
		//set new weight probabilities
		for (k=0; k<K;k++){
			sum_EMG_f 		= getSum(mixture_E[0][k*2], data->N);
			sum_EMG_r 		= getSum(mixture_E[1][k*2], data->N);
			rvs[k].set_weights((sum_EMG_f+sum_EMG_r)/tot );
		}
		t++;
	}
	return 1;
}


int fit(annotations * A,bool DGU, bool EGU, bool verbose, int strand){
	bool FAIL;
	map<string, annotation_cluster *> collections 	= A->collections;
	typedef map<string, annotation_cluster *>::iterator it_type;
	annotation_cluster * cluster;
	annotation * current;
	double threshold 	= 0.0;
	double p 	=0;
	bool curr_strand;
	for(it_type iterator 	= collections.begin(); iterator != collections.end(); iterator++) {
		cluster 			= iterator->second;
		for (cluster = iterator->second; cluster!=NULL; cluster=cluster->next){
			for (current = cluster->root; current!=NULL; current=current->next){
				if (current->density > threshold){
					if (current->strand == "+"){
						curr_strand 	= 0;
					}else if(current->strand=="-"){
						curr_strand 	= 1;
					}

					if (strand == 2 or curr_strand == strand){
						if (DGU){
							FAIL 	= fit_DG(current);
						}else if (EGU){
							FAIL 	= fit_EG(current, 1);
						}
					}
				}
				p+=1;
				if (verbose){
					printf("Running Model Fit: %f percent done\r", (p/A->N)*100);
					cout.flush();
				
				}
				
			}
		} 
		
	}
	cout<<endl;
	return 1;
}
