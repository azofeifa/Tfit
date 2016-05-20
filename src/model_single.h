#ifndef model_single_H
#define model_single_H
#include "load.h"
#include "model.h"

using namespace std;

class NLR{
public:
	double mu, si, l, pi, wn , wl, wr, fp;
	double EX, EX2, EY, EPI, EPIN,   WE, WF, WR;
	double alpha_1, beta_1, alpha_2, beta_2, alpha_3;
	EMG bidir;
	UNI forward;
	UNI reverse;
	NLR();
	void init( segment *, double, int, double,
	 double, double, double, double, double);
	double addSS(double, double, double);
	double pdf(double );
	double get_all();
	void resetSS();
	void set_new_parameters(double);
	void print();

};


class classifier_single{
public:
	double ll, covergence_threshold, max_iterations;
	int K, type;
	double scale;
	NLR * components;
	double alpha_1, beta_1, alpha_2, beta_2, alpha_3;
	double fit(segment *, double);
	classifier_single();
	classifier_single(double, int, int, int, 
		double, double, double, double, double, double);

};


#endif