#ifndef model_H
#define model_H

class UNI{
public:
	double a,b,w,pi;
	
	int st;
	//sufficient stats
	double ri_forward, ri_reverse; //current responsibility
	double r_forward, r_reverse; //running total
	double delta_a, delta_b;
	UNI();
	UNI(double, double, double, int);
	double pdf(double,int);	
	string print();

};
double sum(double * , int );
double LOG(double );


class EMG{
public:
	double mu, si, l, pi, w;
	//sufficient stats
	double ri_forward, ri_reverse; //current responsibility
	double ey, ex, ex2, r_forward, r_reverse;//running total
	EMG();
	EMG(double, double, double, double, double);
	double pdf(double,int);
	double EY(double ,int);
	double EY2(double ,int);
	string print();
};

class NOISE{
public:
	double a, b, w, pi;
	double ri_forward, ri_reverse; //current responsibility
	double r_forward, r_reverse; //running total
	double pdf(double, int);
	NOISE();
	NOISE(double, double, double, double);
};


class component{
public:
	EMG bidir;
	UNI forward;
	UNI reverse;
	NOISE noise; 

	bool type;

	//=====================================
	//parameters to simulate from
	double alpha_0, alpha_1, alpha_2, beta_0, beta_1, beta_2;
	component();
	void initialize(double, double , double, int , double , double, double);

	double evaluate(double, int);
	void add_stats(double, double , int, double);
	void update_parameters(double);
	double get_all_repo();
	void reset();
	void print();
};


class classifier{
public:
	int K; //number of components
	double convergence_threshold; //convergence check
	int max_iterations; //stop after this many iterations
	bool seed; //seed with a gross peak finder
	double noise_max; //fit a uniform noise component, never let it get above this weight
	//===================================================================================
	//Bayesian Priors
	double alpha_0; //symmetric prior for mixing weights; dirchlet
	double beta_0; //symmetric prior for strand probabilities; beta
	double m_0,tau; //priors for component mus; gaussian
	double alpha_1,beta_1; //priors for component sigmas; gamma
	double alpha_2,beta_2; //priors for component lambdas; gamma
	double p;
	int fit(segment *,vector<double>);
	classifier(int);
	component * rvs = NULL;

};

#endif
