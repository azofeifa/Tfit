#include "BIC.h"
#include "model.h"
#include <cmath>
using namespace std;
double BIC3(double ** X, int j, int k, int i,
	    double N_pos, double N_neg,  double sigma, double lambda, double fp, double pi, double w){

  double N                = N_pos + N_neg;
  double l     = X[0][k] - X[0][j];

  double uni_ll= LOG(pi/  (l))*N_pos + LOG((1-pi)/ (l))*N_neg;


  double pi2      = (N_pos+10000) / (N_neg + N_pos+20000);

  double emg_ll   = 0, p1=0.0,p2=0.0;
  EMG EMG_clf(X[0][i], sigma, lambda, w, pi2  );
  EMG_clf.foot_print      = fp;
  
  for (int i = j; i < k;i++ ){
    p1 = EMG_clf.pdf(X[0][i],1) + (1.0-w)*pi*(1.0/l) , p2 = EMG_clf.pdf(X[0][i],-1) + (1.0-w)*(1.0-pi)*(1.0/l) ;
    if (p1 > 0 and p2 > 0 ){//this should always evalulate!!
      emg_ll+=LOG( p1 )*X[1][i] + LOG( p2 )*X[2][i];
    }else{
    }
  }
  double emg_ratio        = (-2*uni_ll + LOG(N)) / (-2*emg_ll + 20*LOG(N))  ;
  //printf("%f,%f,%f, %f, %f,%f,%f\n", sigma, lambda, w, pi2, fp, emg_ratio, l);
  return emg_ratio;
}

