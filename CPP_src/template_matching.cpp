
#include "load.h"
#include "model.h"
#include <math.h> 
#include <limits>
#include <iostream>
#include <algorithm>
#include "template_matching.h"
#include <random>
using namespace std;

double nINF	=-exp(1000);
double INF 	= exp(1000);


//=========================================================
//Peak Finding

vector<vector<double> > find_peaks(vector<double> values, double ** data){
	vector<vector<double> > peaks;
	for (int i = 1; i < values.size()-1; i++ ){
		if (values[i-1] < values[i] && values[i] > values[i+1]){
			vector<double> current(2); 
			current[0] 	= data[0][i];
			current[1] 	= values[i];
			peaks.push_back(current);
		}
	}
	return peaks; 
}

//=========================================================
//BiDirectional Mappings 

double template_ll(segment * data, int center, 
	int start, int stop, double si, double l ){
	double EMG_ll=0;
	double UNI_ll=0;
	double vl 	= 0.5/(data->X[0][stop] - data->X[0][start]);
	EMG EMG_clf(data->X[0][center], si, l, 1.0, 0.5 );
	for (int i = start; i < stop; i++ ){
		UNI_ll+=(LOG(vl)*data->X[1][i] + LOG(vl)*data->X[2][i]);
		EMG_ll+=(LOG(EMG_clf.pdf(data->X[0][i], 1))*data->X[1][i] + LOG(EMG_clf.pdf(data->X[0][i], -1))*data->X[2][i] );
	}
	return UNI_ll/EMG_ll;
	
}
double window_cov(segment * data, int start, int stop){
	double forward 	= 0;
	double reverse 	= 0; 
	for (int i = start; i < stop; i++ ){
		forward+=data->X[1][i];
		reverse+=data->X[2][i];
	}	
	return LOG(forward)+LOG(reverse) ;
}

vector<vector<double>> bubble_sort(vector<vector<double>> X){ //sort vector of vectors by second
	bool changed=true;
	while (changed){
		changed=false;
		for (int i = 0; i < X.size()-1; i++  )	{
			if (X[i][1] < X[i+1][1]){
				vector<double> copy 	= X[i];
				X[i] 					= X[i+1];
				X[i+1] 					= copy;
				changed=true;
			}
		}
	}

	return X;
}

vector<double> window_search(segment * data, double si, double l, bool kind){
	vector<double> values(data->XN);
	int j,k;
	double window = si + (1.0 / l);
	double vl;
	for (int i = 0; i < data->XN; i++){
		j=i;
		while (j < data->XN && (data->X[0][j] - data->X[0][i]) < window){
			j++;
		}
		k=i;
		while (k >0 and (data->X[0][i] - data->X[0][k]) < window){
			k--;
		}
		if (kind){
			vl =	template_ll(data, i, k,j, si, l);
		}else{
			vl =	window_cov(data, k, j);			
		}
		values[i] 	= vl;
	}

	return values;
}
vector<double> take_avg(vector<double> X1, vector<double> X2 ){
	vector<double> avg;
	if (X1.size() != X2.size()){
		cout<<"problem with template matchers...."<<endl;
		return avg;
	}
	for (int i = 0; i < X1.size(); i++){
		avg.push_back((X1[i] + X2[i])/2.);
	}
	return avg;
}


vector<double> peak_bidirs(segment * data){


	vector<double> coverage_values 	= (window_search(data, 1, 0.1, 0));
	vector<double> template_values 	= (window_search(data, 1, 0.1, 1));
	vector<double> average_values 	= take_avg(coverage_values, template_values);
	
	vector<vector<double>> peaks 	= find_peaks(average_values, data->X);
	
	peaks 							= bubble_sort(peaks);
	vector<double> centers ;
	for (int i = 0 ; i < peaks.size(); i++)	{
		centers.push_back(peaks[i][0]);
	}
	return centers;
}

int sample_centers(vector<double> centers, double p){
	random_device rd;
	mt19937 mt(rd());
	default_random_engine generator;
	geometric_distribution<int> distribution(p);
	int i 	= distribution(mt);
	int s 	= centers.size()-1;
	if (i > s){
		return s;
	}
	return i;
}






