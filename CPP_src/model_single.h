#ifndef model_single_H
#define model_single_H
#include "load.h"
using namespace std;
class NORM{
public:
	double mu, si, w;
	double pdf(double);
	NORM();
	NORM(double, double, double);
	string print();
};

class ELON{
public:
	double a,b,w;
	ELON();
	ELON(double, double, double);
	double pdf(double);
	string print();
};

class NLR{
public:
	double mu, si, wn, l, r, wl, wr;
	double EX, EX2, WN, WR, WL;
	NORM loading;
	ELON forward;
	ELON reverse;
	NLR();
	void init(int, int, segment *, double);
	void addSS(double, double, double);
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
	double fit(segment *);
	classifier_single();
	classifier_single(double, int, int, int, double);

};


#endif