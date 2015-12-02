#include "model_selection.h"
#include "template_matching.h"
#include <stdlib.h> 

int BIC(double penality, segment_fits * S){
	double noise_bic, model_bic;
	if (S->M[0]==nINF and S->M[1]==INF){
		return 0;
	}
	if (S->M[1]==nINF){
		S->M[1]=S->M[0];
	}
	noise_bic 	= -2*S->M[0] + log(S->N)*2;
	model_bic 	= -2*S->M[1] + log(S->N)*7*penality;
	if (noise_bic < model_bic){
		return 0;
	}
	return 1;
}

void compute_true_and_false_positive_rates(double penality, vector<segment_fits *> noise_fits,
	vector<segment_fits *> query_fits, double & tp, double & fp ){
	//compute false positives
	for (int i = 0 ; i < noise_fits.size(); i++){
		int s 	=  BIC(penality, noise_fits[i]);
		fp+=s;
	}
	//compute true positives
	for (int i = 0 ; i < query_fits.size(); i++){
		int s 	=  BIC(penality, query_fits[i]);
		tp+=s;
	}
}




double ROC(vector<segment_fits *> noise_fits, vector<segment_fits *> query_fits, 
	double& AUC , double& TP, double& FP, double& TP_at_fp, double& FP_at_tp, double& optimal_penality ){
	double res 	= 500;
	printf("%d\n", noise_fits.size(), query_fits.size());
	double penality_a 	= -1000; //everything is model
	double penality_b 	= 1000; //everthing is noise
	double delta 		= (penality_b - penality_a)/ res;
	double penality;
	double TPN 			= query_fits.size();
	double TNN 			= noise_fits.size();
	double prev_tp = 1.0, prev_fp = 1.0;
	double max_dist = nINF;
	for (int r = 0; r < res; r++){
		penality 		= penality_a + delta*r;
		double tp = 0,fp = 0;
		compute_true_and_false_positive_rates(penality, noise_fits, query_fits, tp, fp);
		tp /= TPN;
		fp /= TNN;
		AUC+=(prev_fp - fp )*prev_tp;
		if (tp -fp   > max_dist   ){
			TP 	= tp;
			FP 	= fp;
			max_dist 	= tp - fp ;
			optimal_penality 	= penality;
		}
		prev_fp=fp, prev_tp=tp;
		if (fp < pow(10,-5)){
			break;
		}
		if (tp > 0.95){
			FP_at_tp=fp;
		}
		if (fp > 0.05){
			TP_at_fp=tp;
		}

	}

}







