#ifndef pdfs_H
#define pdfs_H

class EGU{
public:

	double mu, s, l, w;
	double b, pi_emg, pi_uni;

	
	double pdf_U(double, int);
	double pdf_EG(double, int);
	void set(double, double , double , double , double , double , double);
	void set_partial(double , double , double , double , double );
	void set_weights(double);
	double cond_y(double, int);
	double cond_y_sq(double, int);
	double cond_x(double, int);
	double cond_x_sq(double, int);

	EGU();
};
class DGU{
public:
	double w, i, u,s,l;
	DGU();
	void set(double, double, double, double, double);
	double pdf_U(double);
	double pdf_DG(double);

};
#endif