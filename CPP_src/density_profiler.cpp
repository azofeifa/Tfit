#include "density_profiler.h"
#include <iostream>
#include "split.h"
#include <fstream>
#include <math.h> 
#include <cmath>
#ifdef USING_ICC
#include <aligned_new>
#endif
using namespace std;

gap_interval::gap_interval(){};
gap_interval::gap_interval(double st, double sp){
		start=st, stop=sp;
}

double get_table_mean_var(string GAP_FILE, string bed_file, double res, double bin_res, double scale){


	map<int, vector<double> > A;
	map<string, vector<gap_interval>> G;
	ifstream FH(GAP_FILE);
	if (FH){
		string line, chrom, start, stop;
		vector<string> lineArray;

		while (getline(FH, line)){
			lineArray 	= splitter(line, "\t");
			chrom=lineArray[0], start=lineArray[1], stop=lineArray[2];
			gap_interval I(stof(start), stof(stop));
			G[chrom].push_back(I);
		
		
		}
	}else{
		cout<<"Coudn't open: "<<GAP_FILE<<endl;
	}
	ifstream BED(bed_file);
	if (BED){
		string line, chrom;
		int start, stop;
		double cov;
		vector<string> lineArray;
		string prevchrom="";
		int j,N;
		int t = 0;
		while (getline(BED, line)){
			lineArray 	= splitter(line, "\t");
			chrom=lineArray[0], start=stoi(lineArray[1]), stop=stoi(lineArray[2]), cov=abs(stof(lineArray[3]));
			if (chrom!=prevchrom){
				if (G.find(chrom)!=G.end()){
					j= 0,N=G[chrom].size();
					if (t > 5){
						break;
					}
					t++;

				}else{
					j=0,N=0;
				}

			}
			while (j < N and G[chrom][j].stop<start){
				j++;
			}
			if (j < N and G[chrom][j].start<stop){
				for (int i = start; i < stop; i++){
					G[chrom][j].X.push_back(double(i));
					G[chrom][j].Y.push_back(cov);
					
				}
			}
			prevchrom=chrom;
		}	
	}else{
		cout<<"Coudn't open: "<<bed_file<<endl;
	}

	BED.close();
	//make vector
	typedef map<string, vector<gap_interval> >::iterator it_type;
	map<double, vector <double> > windows;
	for (it_type c = G.begin(); c!=G.end(); c++){
		vector<vector<double>> X;
		for (int g = 0; g <c->second.size(); g++ ){
			for (int i = 0 ; i < c->second[g].X.size(); i++ ){
				vector<double> current(2);
				current[0] 	= c->second[g].X[i], current[1]=c->second[g].Y[i];
				X.push_back(current);
			}
		}
		if (not X.empty()){
			//bin
			double minX = X[0][0];
			double maxX = X[X.size()-1][0];
			int BINS 	 = (maxX-minX)/bin_res;

			double ** XX 	= new double * [BINS];
			for (int j = 0; j < BINS; j++){
				XX[j] 		= new double[2];
				XX[j][0]	= minX + j*bin_res;
				XX[j][1] 	= 0;
			}
			
			int j=0,k=0;
			for (int i = 0; i < X.size(); i++){
				while (j < BINS and XX[j][0] <= X[i][0] ){
					j++;
				}
				if (j < BINS and X[i][0] <=XX[j][0]  ){
					XX[j][1]+=X[i][1];
				}
			}
			double window;
			double window_a = 1000;
			double window_b = 3000;
			double window_d = (window_b-window_a) / res;
			for (int w = 0; w < res; w++){
				window 		= window_a + window_d*w;
				if (windows.find(window) == windows.end()){
					vector<double> curr(2);
					curr[0] 			= 0, curr[1]=0;
					windows[window] 	= curr;
				}
				j=0,k=0;
			
				double x , density ;
				double S 	= 0;
				double AN 	= 0;
				double ALL_SUM=0;
				for (int i = 0 ; i < BINS; i++ ){
					x 	= XX[i][0];
					while (j < BINS and (XX[j][0] - x) < -window  ){
						S-=XX[j][1];
						j++;
					}
					while (0<= k and k < BINS and (XX[k][0]-x) < window){
						S+=XX[k][1];
						k++;
					}
					density 	= S / window;
					windows[window][0]+=density;
					windows[window][1]+=1;
				}
			}
			for (int j = 0 ; j < BINS; j++){
				delete XX[j];
			}
			delete XX;
		}
			
	}

	typedef map<double, vector <double> >::iterator it_type_3;
	double collect = 0;
	double collect_N=0;
	for (it_type_3 w = windows.begin(); w!= windows.end(); w++){
		collect+=(w->second[0]/w->second[1]);
		collect_N++;
	}


	return collect/collect_N;

}


