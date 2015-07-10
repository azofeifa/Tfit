#include <string>
#include <vector>
#include "load.h"
#include "split.h"
#include <iostream>
#include <fstream>
using namespace std;

//========================
//segment class
segment::segment(string chr, int st, int sp){
	chrom	= chr;
	start	= st;
	stop	= sp;
	N 		= 0;
}
segment::segment(){
	N 		= 0;
}


void segment::add(int strand, double x, double y){
	vector<double> v2(2);
	v2[0] 	= x;
	v2[1] 	= y;
	if (forward.empty() && reverse.empty()){
		minX=x;
		maxX=x;
	}else{
		if (x < minX){
			minX=x;
		}
		if (x > maxX){
			maxX=x;
		}
	}
	if (strand == 1){
		forward.push_back(v2);
	}else if (strand==-1){
		reverse.push_back(v2);
	}
}

void segment::bin(int BINS, double scale){
	X 				= new double*[3];
	SCALE 			= scale;
	for (int j = 0 ; j < 3;j++){
		X[j] 		= new double[BINS];
	}
	double delta 	= (maxX - minX) / BINS;
	XN 				= BINS;
	//===================
	//populate bin ranges
	X[0][0] 		= minX;
	for (int i = 1; i < BINS; i++){
		X[0][i] 	= X[0][i-1] + delta;
		X[1][i] 	= 0;
		X[2][i] 	= 0;
	}
	//===================
	//insert forward strand
	int j 	=0;
	for (int i = 0 ; i < forward.size(); i++){
		while (j < BINS and X[0][j] <=forward[i][0]){
			j++;
		}
		X[1][j-1]=X[1][j-1]+forward[i][1];
		N=N+forward[i][1];
	}
	j 	=0;
	//===================
	//insert reverse strand
	for (int i = 0 ; i < reverse.size(); i++){
		while (j < BINS and X[0][j] <=reverse[i][0]){
			j++;
		}
		X[2][j-1]=X[2][j-1]+reverse[i][1];
		N=N+reverse[i][1];
	}
	//===================
	//scale data down for numerical stability
	for (int i = 0; i < BINS; i ++ ){
		X[0][i] 	= (X[0][i]-minX)/scale;
	}
	maxX 			= (maxX-minX)/scale;
	minX 			=0;
	forward.clear();
	reverse.clear();
}

vector<segment*> load_EMGU_format_file(string FILE, string spec){
	ifstream FH(FILE);
	string line; 
	vector<segment*> segments;
	if (FH){
		vector<string> lineArray;
		string chrom, start, stop;
		int strand;
		string x,y;
		segment * S =NULL;
		bool collect=false;
		while (getline(FH, line)){
			if ("#"==line.substr(0,1)){
				if (collect){
					segments.push_back(S);
				}
				lineArray 	= splitter(line.substr(1), ",");
				chrom 		= lineArray[0];
				collect		= false;
				if (spec=="all" or chrom==spec){
					start 		= lineArray[1];
					stop 		= lineArray[2];
					S 			= new segment(chrom, stoi(start), stoi(stop) );
					collect 	= true;
				}
			}else if ("~"==line.substr(0,1) and collect){
				if (line.substr(1)=="forward"){
					strand 		= 1;
				}else{
					strand 		= -1;
				}
			}else if (collect){
				lineArray 	= splitter(line, ",");
				x 			= lineArray[0];
				y 			= lineArray[1];
				S->add(strand, stod(x), stod(y) );
			}
		}
	if (collect){
		segments.push_back(S);
	}
	
	}else{
		cout<<"could not open: "<<FILE<<" "<<endl;
		return segments;
	}
	return segments;
}


void BIN(vector<segment*> segments, int BINS, double scale){
	for (int i = 0 ; i < segments.size() ; i ++){
		segments[i]->bin(BINS, scale);
	}
}
