#include <string>
#include <vector>
#include "load.h"
#include "split.h"
#include <iostream>
#include <fstream>
#include <map>
#include "dirent.h"
#include "template_matching.h"
#include "model.h"
#include "read_in_parameters.h"
#include "across_segments.h"
#include "model_selection.h"
#include <cmath>
#include <math.h> 
#include <limits>
#include <iostream>
#include <algorithm>
#include <unistd.h>
#include <random>

#include <stdio.h>
#include <time.h>
#include <math.h> 

using namespace std;
//========================
//model fit wrapper classes
model_component::model_component(){}

model_component::model_component(simple_c sc){
	mu 	= sc.ps[0],si = sc.ps[1], l = sc.ps[2], w_e = sc.ps[3], pi = sc.ps[4];
	f_b = sc.ps[7], f_w = sc.ps[5], r_a = sc.ps[8], r_w = sc.ps[6];
}

bidir_preds::bidir_preds(double ll, double NN){
	noise_ll 	= ll;
	N=NN;
}

void bidir_preds::model_selection(double penality){
	typedef map<int, all_model_components >::iterator it_type;
	BIC_score 	= -2*noise_ll + 1*log(N);
	double curr_BIC_score;
	argK 	= 0;
	printf("--------------------------\n");
	printf("Noise LL: %f, Noise BIC Score: %f, N: %f\n", noise_ll, BIC_score, N );
	for (it_type s = G.begin(); s!=G.end(); s++){
		curr_BIC_score 	= -2*s->second.ll + penality*s->first*7*log(N);
		if (curr_BIC_score < BIC_score){
			BIC_score 	= curr_BIC_score, argK 	= s->first;
		}
		printf("K: %d, LL: %f, Noise BIC Score: %f\n", s->first
		, s->second.ll, curr_BIC_score  );
	

	}
	if (argK!= 0){
		arg_all_model_components 	= G[argK];
	
	}
	printf("BEST: %d\n", argK);
	printf("--------------------------\n");
}

bidir_preds::bidir_preds( ){
	noise_ll=0, N=0;
}

void bidir_preds::insert_component(int K, simple_c sc, double LL){
	G[K].insert_component(sc, LL);
}

all_model_components::all_model_components(){
	ll 	= nINF;
}

void all_model_components::insert_component(simple_c sc, double LL){
	ll=LL;
	all_components.push_back(model_component(sc));
}

segment::segment(string chr, int st, int sp){
	chrom	= chr;
	start	= st;
	stop	= sp;
	N 		= 0;
	minX=st, maxX=sp;
	counts 	= 1;
	XN 		= 0;
	ID 		= 0;
	strand 	= ".";
	chrom_ID= 0;

}
segment::segment(string chr, int st, int sp, int i){
	chrom	= chr;
	start	= st;
	stop	= sp;
	N 		= 0;
	minX=st, maxX=sp;
	counts 	= 1;
	XN 		= 0;
	ID 		= i;
	strand 	= ".";
	chrom_ID= 0;

}

segment::segment(string chr, int st, int sp, int i, string STR){
	chrom	= chr;
	start	= st;
	stop	= sp;
	N 		= 0;
	minX=st, maxX=sp;
	counts 	= 1;
	XN 		= 0;
	ID 		= i;
	strand 	= STR;
	chrom_ID= 0;

}

segment::segment(){
	N 		= 0;
	counts 	= 1;
	XN 		= 0;
	ID 		= 0;
	strand 	= ".";
	chrom_ID= 0;
}

string segment::write_out(){
	string text 	= ("#" + chrom + ":" + to_string(start) + "-" 
		+ to_string(stop) + "," + to_string(int(N))+ "\n");
	return text;
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
	N+=(x+y);
}

void segment::add2(int strand, double x, double y){
	vector<double> v2(2);
	v2[0] 	= x;
	v2[1] 	= y;
	if (forward.empty() && reverse.empty()){
		minX=x;
		maxX=x;
	}else{
		if (x < minX){
			minX=x;
			start=int(x);
		}
		if (x > maxX){
			maxX=x;
			stop=int(x);
		}
	}
	if (strand == 1){
		forward.push_back(v2);
	}else if (strand==-1){
		reverse.push_back(v2);
	}
}

void segment::add_fitted_bidir(vector<double> D){
	fitted_bidirs.push_back(D);
	start 	= min(start, int(D[0]));
	stop 	= max(stop, int(D[1]));
}


void segment::bin(double delta, double scale, bool erase){
	X 				= new double*[3];
	SCALE 			= scale;
	int BINS;
	BINS 		= (maxX-minX)/delta;
	start = minX, stop=maxX;
	for (int j = 0 ; j < 3;j++){
		X[j] 		= new double[BINS];
	}
	N 				= 0;
	XN 				= BINS;
	//===================
	//populate bin ranges
	X[0][0] 		= double(minX);
	X[1][0]=0,X[2][0]=0;
	for (int i = 1; i < BINS; i++){
		X[0][i] 	= X[0][i-1] + delta;
		X[1][i] 	= 0;
		X[2][i] 	= 0;
	}
	// ===================
	//insert forward strand
	int j 	=0;
	//printf("start: %d , stop: %d , bins: %d ,delta: %f, forward: %d, reverse: %d\n", start, stop, BINS, delta, forward.size(), reverse.size() );
	for (int i = 0 ; i < forward.size(); i++){
		while (j < BINS and X[0][j] <=forward[i][0]){
			j++;
		}
		if (j < BINS and forward[i][0]<= X[0][j]){
			X[1][j-1]+=forward[i][1];
			N+=forward[i][1];
		}
	}
	j 	=0;
	//===================
	//insert reverse strand
	for (int i = 0 ; i < reverse.size(); i++){
		while (j < BINS and X[0][j] <=reverse[i][0]){
			j++;
		}
		if (j < BINS and reverse[i][0]<= X[0][j]){
			X[2][j-1]+=reverse[i][1];
			N+=reverse[i][1];
		}
	}
	//===================
	//scale data down for numerical stability
	if (scale){
		for (int i = 0; i < BINS; i ++ ){

			X[0][i] 	= (X[0][i]-minX)/scale;
		}
	}
	//we also want to get rid of those data points that we don't need
	//i.e. the ones where there is no data coverage values on either the 
	//forward or reverse strands

	int realN 		= 0;
	for (int i = 0; i < BINS;i++){
		if (X[1][i]>0 or X[2][i]>0){
			realN++;
		}
	}
	if (erase){
		double ** newX 	= new double*[3];
		for (int j=0; j<3;j++){
			newX[j] 	= new double[realN];
		}
		j = 0;
		for (int i = 0; i < BINS; i ++){
			if (X[1][i]>0 or X[2][i]>0){
				newX[0][j] 	= X[0][i];
				newX[1][j] 	= X[1][i];
				newX[2][j] 	= X[2][i];
				j++;
			}
		}
		if (realN!=j){
			printf("WHAT? %d,%d\n", j, realN);
		}
		//clear previous memory
		for (int i = 0; i < 3; i ++){
			delete X[i];
		}
		delete X;
		X 				= newX;
		XN 				= realN;
	}
	if (scale){
		if (not centers.empty()){
			for (int i = 0; i < centers.size(); i++){
				centers[i]=(centers[i]-minX)/scale;			
			}
		}
		if (not fitted_bidirs.empty() ){
			for (int fb = 0; fb < fitted_bidirs.size(); fb++){
				int center 	= fitted_bidirs[fb][0];
				int std 	= fitted_bidirs[fb][1]*0.5 + (1. /  fitted_bidirs[fb][2]);
				int a 		= center - std*3;
				int b 		= center + std*3;
				


				fitted_bidirs[fb][0] = (fitted_bidirs[fb][0] - minX)/scale;
				fitted_bidirs[fb][1] /= scale;
				fitted_bidirs[fb][2] *= scale; 
			}
		}

		maxX 			= (maxX-minX)/scale;
		minX 			=0;
	}
	double S=0;
	for (int i = 0; i < XN; i++){
		S+=X[1][i];
	}
	forward.clear();
	reverse.clear();
}

void segment::insert_bidirectional_data(int pad){

	for (int i = 0; i < bidirectional_bounds.size();i++){
		bidirectional_bounds[i][0]= max(bidirectional_bounds[i][0]-pad, minX);
		bidirectional_bounds[i][1]= min(bidirectional_bounds[i][1]+pad, maxX);
	}
	vector<vector<double>> merged_bidirectional_bounds;

	//first thing merge
	int bound_N 		= bidirectional_bounds.size();
	int i 	= 0;
	double o_st, o_sp;
	int cts  	= 0;
	while (i < bound_N){
		o_st  	= bidirectional_bounds[i][0], o_sp=bidirectional_bounds[i][1];
		cts 	= 0;
		while (i < bound_N and 
			( bidirectional_bounds[i][1] > o_st) and 
			( bidirectional_bounds[i][0] < o_sp)  ){

			o_st=min(o_st, bidirectional_bounds[i][0]);
			o_sp=min(o_sp, bidirectional_bounds[i][1]);
			i++;
			cts++;
		} 
		vector<double> merged(2);
		merged[0]=o_st, merged[1]=o_sp;
		merged_bidirectional_bounds.push_back(merged);
		bidir_counts.push_back(cts);
	
	}
	//first thing is to get NS
	vector<double> Y_f(2);
	vector<double> Y_r(2);
	
	segment * bidir_seg 	= NULL;
	for (int j = 0; j < merged_bidirectional_bounds.size(); j++ ){
		bidir_seg 	= new segment(chrom,
				merged_bidirectional_bounds[j][0],
				merged_bidirectional_bounds[j][1] );
		bidir_seg->counts 	= bidir_counts[j];
		bidirectional_data.push_back(bidir_seg);
	}



	int j 	= 0;
	for (int i = 0; i < XN;i++){
		while (j < bidirectional_data.size() and
			bidirectional_data[j]->stop < X[0][i] ){
			j++;
		
		}
		if (j < bidirectional_data.size() and bidirectional_data[j]->start <=X[0][i] and  
			X[0][i]<=bidirectional_data[j]->stop){
		
				Y_f[0]=X[0][i],Y_f[1]=X[1][i];
				Y_r[0]=X[0][i],Y_r[1]=X[2][i];
				bidirectional_data[j]->forward.push_back(Y_f);
				bidirectional_data[j]->reverse.push_back(Y_r);
		}
	}
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

vector<vector<double> > interval_sort(vector<vector<double> > A){
	bool GOOD=false;
	vector<double> cp;
	if (not A.size()){
		return A;
	}
	while (not GOOD){
		GOOD=true;
		for (int i =1; i < A.size(); i++){
			if (A[i][0] < A[i-1][0]){
				cp 		= A[i-1];
				A[i-1] 	= A[i];
				A[i] 	= cp;
				GOOD=false;
			}
		}
	}
	return A;
}

vector<segment *> bubble_sort_segment(vector<segment *> A){
	bool GOOD=false;
	segment * cp;
	if (not A.size()){
		return A;
	}
	while (not GOOD){
		GOOD=true;
		for (int i =1; i < A.size(); i++){
			if (A[i]->start < A[i-1]->start){
				cp 		= A[i-1];
				A[i-1] 	= A[i];
				A[i] 	= cp;
				GOOD=false;
			}
		}
	}
	return A;	
}


void BIN(vector<segment*> segments, int BINS, double scale, bool erase){
	for (int i = 0 ; i < segments.size() ; i ++){
		if (segments[i]->forward.size() > 0 or segments[i]->reverse.size() > 0 ){
			segments[i]->bin(BINS, scale, erase);
		}
	}
}

interval::interval(){
	hits 	= 0;
	EMPTY 	= true;
};

interval::interval(string chr, int st, int sp){
	chrom 	= chr;
	start 	= st;
	stop 	= sp;
	hits 	= 0;
	EMPTY 	= false;
};
interval::interval(string chr, int st, int sp, int IDD){
	chrom 	= chr, start =st , stop = sp, ID = IDD, EMPTY=false;
}
interval::interval(string chr, int st, int sp, int IDD, string str, vector<vector<double>> p ){
	
	chrom 	= chr, start =st , stop = sp, ID = IDD, EMPTY=false, STRAND=str;
	parameters 	= p;
}
interval::interval(string chr, int st, int sp, int IDD, string str, vector<vector<double>> p, int ct ){
	
	chrom 	= chr, start =st , stop = sp, ID = IDD, EMPTY=false, STRAND=str;
	parameters 	= p;
	counts 	= ct;
}


void interval::insert(double x, double y, int strand){
}

vector<interval> bubble_sort(vector<interval> X){
	bool changed=true;
	while (changed){
		changed=false;
		for (int i = 0; i < X.size()-1; i++  )	{
			if (X[i].start > X[i+1].start){ //sort by starting position
				interval copy 			= X[i];
				X[i] 					= X[i+1];
				X[i+1] 					= copy;
				changed=true;
			}
		}
	}
	return X;
}


merged_interval::merged_interval(){};

merged_interval::merged_interval(int st, int sp, interval unit, int i){
	start 	= st;
	stop 	= sp;
	intervals.push_back(unit);
	id 		= i;
};

void merged_interval::update(interval insert){
	if (insert.start!= start and insert.stop != stop){
		start= min(start, insert.start), stop= max(stop, insert.stop);
		intervals.push_back(insert);
	}
}

void merged_interval::insert(double x,double y, int strand){
	if (strand == 1){
		for (int i=0; i < intervals.size(); i++){
			if (intervals[i].start <= x and x<= intervals[i].stop){
				intervals[i].forward_x.push_back(x);
				intervals[i].forward_y.push_back(y);
			}
		}
	}else{
		for (int i=0; i < intervals.size(); i++){
			if (intervals[i].start <= x and x<= intervals[i].stop){
				intervals[i].reverse_x.push_back(x);
				intervals[i].reverse_y.push_back(y);
			}
		}
	}
}

bool merged_interval::find(int ST, int SP){
	bool FOUND=false;
	for (int i = 0; i < intervals.size(); i++){
		if (  SP > intervals[i].start and   ST < intervals[i].stop){
			intervals[i].hits+=1;
			FOUND=true;
		}
	}
	return FOUND;
}
interval merged_interval::get_interval(int ST, int SP){
	bool FOUND=false;
	interval empty;
	for (int i = 0; i < intervals.size(); i++){
		if (  ST == intervals[i].start and   SP == intervals[i].stop){
			return intervals[i];
		}
	}
	return empty;
}

int merged_interval::get_hits(bool y, int ST, int SP){
	int HITS = 0;
	for (int i = 0; i < intervals.size(); i++){
		if ( SP > intervals[i].start  and ST< intervals[i].stop ){
			if (not y ){
				HITS+=intervals[i].hits;
			}else if (intervals[i].hits>0){
				HITS+=1;
			}
		}
	}
	return HITS;
}

void merged_interval::reset_hits(){
	for (int i = 0; i < intervals.size(); i++){
		intervals[i].hits=0;
	}
}

int merged_interval::get_total(int ST, int SP){
	if (SP > start and ST < stop){
		return intervals.size();
	}
	return 0.;
}



//================================================================================================
//interval tree code

interval_tree::interval_tree(){
	left 	= NULL;
	right 	= NULL;
	current = NULL;
};

int interval_tree::getDepth(){
	if (left!=NULL and right!=NULL){
		return 1 + left->getDepth() + right->getDepth();
	}else if(left!=NULL){
		return 1 +  left->getDepth();
	}else if (right!=NULL){
		return 1 +  right->getDepth();	
	}
	return 1;
}

int interval_tree::get_total(int start, int stop){

	if (left!=NULL and right!=NULL){
		return current->get_total(start, stop) + left->get_total(start, stop) + right->get_total(start, stop);
	}else if(left!=NULL){
		return current->get_total(start, stop) +  left->get_total(start, stop);
	}else if (right!=NULL){
		return current->get_total(start, stop) +  right->get_total(start, stop);	
	}
	return current->get_total(start, stop);
}

void interval_tree::construct(vector<merged_interval *> D){
	if (D.size()>0){
		int N 		= D.size();
		int center 	= N / 2;
		current 	= D[center];
		ID 			= D[center]->id;
		//printf("%d,%d, %d,%d\n", ID, current->start, N, center);
		if ((center ) > 0){
			left 	= new interval_tree();
			vector<merged_interval *> LFT(D.begin(), D.begin() + center);
			left->construct(LFT);
		}else{
			left 	= NULL;
		}
		if (center + 1 < N){
			right 	= new interval_tree();
			vector<merged_interval *> RT(D.begin()+center, D.end()  );
			right->construct(RT);
		}else{
			right 	= NULL;
		}
	}
}
void interval_tree::insert(double x, double y, int strand){
	if (x < current->start and left!=NULL){
		left->insert(x,y,strand);
	}else if (x > current->stop and right !=NULL){
		right->insert(x,y,strand);
	}else if (current!=NULL){
		current->insert(x,y,strand);
	}
}
void interval_tree::insert_into_array(merged_interval * ARRAY , int N ){
	ARRAY[ID] 	= *current;
	if (left!=NULL){
		left->insert_into_array(ARRAY, N);
	}
	if (right != NULL){
		right->insert_into_array(ARRAY, N);
	}
}

int interval_tree::get_hits(bool y, int start, int stop){
	if (left!=NULL and right!=NULL){

		return current->get_hits(y,start, stop) + left->get_hits(y, start, stop) + right->get_hits(y,start, stop);
	}else if(left!=NULL){
		return current->get_hits(y,start, stop) +  left->get_hits(y,start, stop);
	}else if (right!=NULL){
		return current->get_hits(y,start, stop) +  right->get_hits(y,start, stop);	
	}
	return 0;
}

void interval_tree::reset_hits( ){
	current->reset_hits();
	if (left!=NULL and right!=NULL){
		left->reset_hits();
		right->reset_hits();
	}else if(left!=NULL){

		left->reset_hits();
	}else if (right!=NULL){
		right->reset_hits();	
	}
}

bool interval_tree::find(int start, int stop){
	if (( stop > current->start) and ( start < current->stop) ){
		return current->find(start, stop);
	}
	else if (stop < current->start and left!=NULL){
		return left->find(start, stop);
	}
	else if (start > current->stop and right!=NULL){
		return right->find(start, stop);
	}else{
		return false;
	}
}


interval interval_tree::get_interval(int start, int stop){
	interval empty;
	if (( stop > current->start) and ( start < current->stop) ){
		return current->get_interval(start, stop);
	}
	else if (stop <= current->start and left!=NULL){
		return left->get_interval(start, stop);
	}
	else if (start >= current->stop and right!=NULL){
		return right->get_interval(start, stop);
	}
	return empty;

}

//================================================================================================
//loading from file functions...need to clean this up...

vector<segment*> load_bedgraphs_total(string forward_strand, 
	string reverse_strand, int BINS, double scale, string spec_chrom, map<string, int>& chromosomes
	, map<int, string>& ID_to_chrom){
	bool FOUND 	= false;
	if (spec_chrom=="all"){
		FOUND 	= true;
	}
	map<string, segment*> 	G;
	vector<segment*> segments;
	
	ifstream FH(forward_strand);
	ifstream FH2(reverse_strand);
	if (not FH){
		cout<<"couldn't open "<<forward_strand<<endl;
		return segments;
	}else if (not FH2){
		cout<<"couldn't open "<<reverse_strand<<endl;
		return segments;
	}
	
	string line, chrom;
	int start, stop;
	double coverage;
	vector<string> lineArray;
	string prevChrom="";
	segment * S =NULL;
	bool INSERT 	= false;
	while (getline(FH, line)){
		lineArray=splitter(line, "\t");
		chrom=lineArray[0], start=stoi(lineArray[1]), stop=stoi(lineArray[2]), coverage=abs(stof(lineArray[3]));
		if (chrom != prevChrom and (chrom==spec_chrom or spec_chrom=="all")  )  {
			FOUND 		= true;
			if (chrom.size() < 6){
				G[chrom] 	= new segment(chrom, start, stop );
				INSERT 		= true;
			}else{
				INSERT 		= false;
			}
		}
		if (FOUND and chrom!= spec_chrom and spec_chrom!= "all"){
			break;
		}
		if (INSERT){
			G[chrom]->add2(1, double((stop + start) / 2.), coverage);
		}
		prevChrom=chrom;

	}


	prevChrom="";
	FOUND=false;
	int c =1;
	while (getline(FH2, line)){
		lineArray=splitter(line, "\t");
		chrom=lineArray[0], start=stoi(lineArray[1]), stop=stoi(lineArray[2]), coverage=abs(stof(lineArray[3]));

		if (FOUND and chrom!= spec_chrom and spec_chrom!= "all"){
			break;
		}
		if (G.find(chrom)!= G.end()){
			G[chrom]->add2(-1,double((stop + start) / 2.), coverage);
			FOUND= true;
		}
		prevChrom=chrom;
	}
	typedef map<string, segment*>::iterator it_type;
	for (it_type i = G.begin(); i != G.end(); i++){
		i->second->bin(BINS, scale,false);
		if (chromosomes.find(i->second->chrom)==chromosomes.end()){
			chromosomes[i->second->chrom]=c;
			ID_to_chrom[c] 	= i->second->chrom;
			c++;
		}


		segments.push_back(i->second);
	}
	if (not FOUND){
		segments.clear();
		printf("couldn't find chromosome %s in bedgraph files\n", spec_chrom.c_str());
	}
	return segments;
}

vector<segment*> load_bedgraphs_single(string forward_strand,
 int BINS, double scale, string spec_chrom){
	bool FOUND 	= false;
	if (spec_chrom=="all"){
		FOUND 	= true;
	}
	map<string, segment*> 	G;
	vector<segment*> segments;
	
	ifstream FH(forward_strand);
	if (not FH){
		cout<<"couldn't open "<<forward_strand<<endl;
		return segments;
	}
	
	string line, chrom;
	int start, stop;
	double coverage;
	vector<string> lineArray;
	string prevChrom="";
	segment * S =NULL;
	
	while (getline(FH, line)){
		lineArray=splitter(line, "\t");
		chrom=lineArray[0], start=stoi(lineArray[1]), stop=stoi(lineArray[2]), coverage=stod(lineArray[3]);
		if (chrom != prevChrom){
			if (chrom == spec_chrom){
				FOUND 	= true;
			}
			G[chrom] 	= new segment(chrom, start, stop );
		}
		G[chrom]->add2(1, double((stop + start) / 2.), coverage);
		prevChrom=chrom;

	}

	typedef map<string, segment*>::iterator it_type;
	for (it_type i = G.begin(); i != G.end(); i++){
		i->second->bin(BINS, scale,false);
		segments.push_back(i->second);
	}
	if (not FOUND){
		segments.clear();
		printf("couldn't find chromosome %s in bedgraph files\n", spec_chrom.c_str());
	}
	return segments;
}

map<string, vector<merged_interval*> > segments_to_merged_intervals(map<string, vector<segment *> > FSI, ofstream& FHW){
	map<string, vector<interval>> G;	
	typedef map<string, vector<segment *> >::iterator it_type_2;
	for (it_type_2 i = FSI.begin(); i != FSI.end(); i++){
		for (int j = 0; j < i->second.size(); j++){
			G[i->first].push_back(interval(i->second[j]->chrom, 
				i->second[j]->start, i->second[j]->stop, 
				i->second[j]->ID, i->second[j]->strand, i->second[j]->parameters, i->second[j]->counts)  );
		}
	}
	//want to sort intervals by there ending point
	map<string, vector<merged_interval*> > A;
	typedef map<std::string, vector<interval>>::iterator it_type;
	int N, i, j;

	for(it_type c = G.begin(); c != G.end(); c++) {
		if (G[c->first].size()>0){
		    G[c->first] 	= bubble_sort(G[c->first]);

		    j 				= 0;
		    merged_interval * I = new  merged_interval(G[c->first][0].start, G[c->first][0].stop, G[c->first][0], j);
		    I->parameters 	= G[c->first][0].parameters;
		    N 				= G[c->first].size();
		    i 				= 1;
		    bool inserted 	= false;

		    while (i < N and not inserted){
		    	while (i<N and G[c->first][i].start <= I->stop and G[c->first][i].stop >= I->start ){
		    		I->update(G[c->first][i]);
		    		i++;
		    	}
		    	A[c->first].push_back(I);
		    	inserted=true;
		    	if (i < N){
		    		I 	= new merged_interval(G[c->first][i].start, G[c->first][i].stop, G[c->first][i], j);
				    I->parameters 	= G[c->first][i].parameters;

		    		inserted=false;
		    		j++;
		    	}
		    //	i++;
		    }
			if (not inserted){
				A[c->first].push_back(I);		    	
			}
		}
	}
	return A;
	
}

map<string, vector<merged_interval*> >  load_intervals(string FILE, int pad){

	bool header = 1;
	ifstream FH(FILE);
	map<string, vector<interval>> G;	
	if (FH){
		string line, chrom, start, stop;
		vector<string> lineArray;

		while (getline(FH, line)){
			if (not header){
				lineArray 	= splitter(line, "\t");
				chrom=lineArray[0], start=lineArray[1], stop=lineArray[2];
				G[chrom].push_back(interval(chrom, stoi(start)-pad, stoi(stop)+pad ));
			
			}else{
				header 		= 0;
			}
		}
	}else{
		cout<<"Coudn't open: "<<FILE<<endl;

	}
	//want to sort intervals by there ending point
	map<string, vector<merged_interval*> > A;
	typedef map<std::string, vector<interval>>::iterator it_type;
	int N, i, j;

	for(it_type c = G.begin(); c != G.end(); c++) {
		if (G[c->first].size()>0){
		    G[c->first] 	= bubble_sort(G[c->first]);

		    j 				= 0;
		    merged_interval * I = new  merged_interval(G[c->first][0].start, G[c->first][0].stop, G[c->first][0], j);
		    N 				= G[c->first].size();
		    i 				= 1;
		    bool inserted 	= false;
		    while (i < N){
		    	while (i<N and G[c->first][i].start < I->stop and G[c->first][i].stop > I->start ){
		    		I->update(G[c->first][i]);
		    		i++;
		    	}
		    	A[c->first].push_back(I);
		    	inserted=true;
		    	if (i < N){
		    		I 	= new merged_interval(G[c->first][i].start, G[c->first][i].stop, G[c->first][i], j);
		    		inserted=false;
		    		j++;
		    	}
		    	i++;
		    }
			if (not inserted){
				A[c->first].push_back(I);
			}
		}
	}



	return A;
}

void insert_bedgraph(map<string, interval_tree *> intervals, string FILE, int strand){

	ifstream FH(FILE);
	if (FH){
		string line, chrom,  prevChrom;
		prevChrom="";
		int j = 0;
		int N = 0;
		int k;
		vector<string> lineArray;
		int start, stop;
		double coverage;
		while (getline(FH, line)){
			lineArray=splitter(line, "\t");
			chrom=lineArray[0], start=stoi(lineArray[1]), stop=stoi(lineArray[2]), coverage=stod(lineArray[3]);
			for (int i=start; i < stop;i++){
				if (intervals.find(chrom) !=intervals.end()){
					
					intervals[chrom]->insert(double(i), coverage, strand);
				}
			}		
		}
	}else{
		cout<<"Coudn't open: "<<FILE<<endl;
	}
}

map<string, vector<interval> > insert_bed_file(string FILE, 
	map<string, vector< interval >> G, string spec_chrom){
	fstream FH(FILE);
	if (FH){
		string line, chrom;
		int start, stop ;
		vector<string> lineArray;
		merged_interval * mi 	= NULL;
		while (getline(FH, line)){
			lineArray 	= splitter(line, "\t");
			start=stoi(lineArray[1]), stop=stoi(lineArray[2]);
			chrom 		= lineArray[0];
			if (chrom==spec_chrom){
				G[chrom].push_back(interval(chrom,start, stop));
			}
		}
	}
	return G;
}

map<string, interval_tree *> load_bidir_bed_files(string in_directory, string spec_chrom){
	DIR *dir;
	struct dirent *ent;
	string FILE;
	map<string, vector<interval>> G;
	map<string, interval_tree *> I;
	if ((dir = opendir (in_directory.c_str())) != NULL) {
		/* print all the files and directories within directory */
		while ((ent = readdir (dir)) != NULL) {
			FILE 	= string(ent->d_name);
			G 		= insert_bed_file(in_directory+FILE, G, spec_chrom);	
		}
	}else {
		/* could not open directory */
		cout<<"couldn't open directory: "<<in_directory<<endl;
		return I;
	}
	if (G.empty()){
		cout<<"couldn't "<<spec_chrom<<"in bed files in directory"<<endl;
		return I;
	}


	closedir (dir);
	//want to sort intervals by there ending point
	map<string, vector<merged_interval*> > A;
	typedef map<std::string, vector<interval>>::iterator it_type;
	typedef map<std::string, vector<merged_interval*>>::iterator it_type_2;
	int N, i, j;
	int u 	= 0;
	for(it_type c = G.begin(); c != G.end(); c++) {
		if (G[c->first].size()>0){
			
		    G[c->first] 	= bubble_sort(G[c->first]);

		    j 				= 0;
		    merged_interval * I = new  merged_interval(G[c->first][0].start, G[c->first][0].stop, G[c->first][0], j);
		    N 				= G[c->first].size();
		    i 				= 1;
		    bool inserted 	= false;
		    while (i < N){
		    	while (i<N and G[c->first][i].start < I->stop and G[c->first][i].stop > I->start ){
		    		I->update(G[c->first][i]);
		    		i++;
		    	}
		    	A[c->first].push_back(I);
		    	inserted=true;
		    	if (i < N){
		    		I 	= new merged_interval(G[c->first][i].start, G[c->first][i].stop, G[c->first][i], j);
		    		inserted=false;
		    		j++;
		    	}
		    	i++;
		    }
			if (not inserted){
				A[c->first].push_back(I);
			}
		}
	}
	for (it_type_2 c = A.begin(); c!=A.end(); c++ ){
		I[c->first] 	= new interval_tree;
		I[c->first]->construct(c->second);
	}



	return I;
}

vector<segment *> bidir_to_segment(map<string , vector<vector<double> > > G, 
	string forward_file, string reverse_file, int pad, string spec_chrom, map<string, int> chromosomes){
	
	typedef map<string , vector<vector<double> > >::iterator it_type;
	typedef map<string, vector<segment *> >::iterator it_type_2;
	
	vector<segment*> segments;
	map<string, vector<segment*> > A;
	map<string, vector<segment*> > B;
	segment * S;
	for (it_type i = G.begin(); i != G.end(); i++){
		G[i->first] 	= interval_sort(i->second);
		for (int j = 0; j < G[i->first].size(); j++){
			S 			= new segment(i->first, int(G[i->first][j][0] )-pad, int(G[i->first][j][1])+pad);
			S->chrom_ID = chromosomes[S->chrom];
			B[i->first].push_back(S);
		}
	
	}

	int N,j;
	//we want to merge the overlaping calls
	int overlaps 	= 0;
	for (it_type_2 i = B.begin(); i!=B.end();i++){
		N 	= i->second.size(), j 	= 1;
		S 	= i->second[0];
		S->centers.push_back((S->start + S->stop) /2.);
		if (N==1){
			A[i->first].push_back(S);				
		}
		while (j < N){
			while (j < N and i->second[j]->start < S->stop and i->second[j]->stop > S->start  ){
				double center1=(S->start + S->stop) /2.;
				double center2=(i->second[j]->start + i->second[j]->stop) /2.;
				
				double dist = abs(center2-center1);
				S->centers.push_back(center2);
				if (dist > 1000){
					S->counts+=1;
				}
				S->start 	= min(i->second[j]->start, S->start), S->stop 	= max(i->second[j]->stop,S->stop );
				j++;
				overlaps++;
			}

			A[i->first].push_back(S);
			if (j < N){
				S 			= i->second[j];
			}
		}
	
	}
	vector<string> FILES = {forward_file, reverse_file};
	int strands[2] 		= {1,-1};
	int start, stop;
	double coverage;
	N 	= 0,j 	= 0;
	int strand 	= 1;
	int o_st, o_sp;
	vector<string> lineArray;
	string chrom, prevchrom, line;
	bool FOUND;
	bool ALL 	= spec_chrom == "all";
	for (int i =0;i<2;i++){
		strand 	= strands[i];
		ifstream FH(FILES[i]);
		if (FH){
			FOUND 	= false;
			prevchrom="";
			while (getline(FH, line)){
				lineArray 	= splitter(line, "\t");
				chrom 		= lineArray[0];
				if (FOUND and not ALL and chrom != spec_chrom  ){
					break;

				}
				start=stoi(lineArray[1]),stop=stoi(lineArray[2]), coverage = abs(stod(lineArray[3]));
				if (prevchrom!=chrom){
					if (A.find(chrom)!=A.end()){
						FOUND=true;
						N 	= A[chrom].size(), j = 0;
					}else{
						N 	= 0,j=0;
					}
				}
				while (j < N and A[chrom][j]->stop < start){
					j++;
				}
				if (j < N and A[chrom][j]->start < stop){ //overlap!
					o_st 	= max(A[chrom][j]->start, start);
					o_sp 	= min(A[chrom][j]->stop, stop);
					if (o_st == o_sp){
						o_sp+=1;
					}
					for (int u = o_st; u < o_sp; u++){
						A[chrom][j]->add(strand, double(u), coverage );
					}
				}
				prevchrom 	= chrom;
			}
			
		}else{
			cout<<"could not open: "<<FILES[i]<<endl;
			segments.clear();

			return segments;
		}
		FH.close();
	}
	for (it_type_2 i = A.begin(); i!=A.end();i++){
		for (int j = 0; j < i->second.size();j++){
			segments.push_back(i->second[j]);
		}
	}
	return segments;
}

vector<segment *> insert_bedgraph_to_segment(map<string, vector<segment *> > A, string forward_file, string reverse_file, int rank){
	vector<string> FILES = {forward_file, reverse_file};
	int strands[2] 		= {1,-1};
	int start, stop, N, j;
	double coverage;
	N 	= 0,j 	= 0;
	int strand 	= 1;
	int o_st, o_sp;
	vector<string> lineArray;
	string chrom, prevchrom, line;
	vector<segment *> segments;
	for (int i =0;i<2;i++){
		strand 	= strands[i];
		ifstream FH(FILES[i]);
		if (FH){
			prevchrom="";
			while (getline(FH, line)){
				lineArray 	= splitter(line, "\t");
				chrom 		= lineArray[0];
				start=stoi(lineArray[1]),stop=stoi(lineArray[2]), coverage = abs(stod(lineArray[3]));

				if (prevchrom!=chrom){
					if (A.find(chrom)!=A.end()){
						N 	= A[chrom].size(), j = 0;
					}else{
						N 	= 0,j=0;
					}
				}
				while (j < N and A[chrom][j]->stop < start){
					j++;
				}
				if (j < N and A[chrom][j]->start < stop){ //overlap!
					o_st 	= max(A[chrom][j]->start, start);
					o_sp 	= min(A[chrom][j]->stop, stop);
					for (int u = o_st; u < o_sp; u++){
						A[chrom][j]->add(strand, double(u), coverage );
					}
				
				}
				prevchrom 	= chrom;
			}
			FH.close();
			
		}else{
			cout<<"could not open: "<<FILES[i]<<endl;
			segments.clear();
			return segments;
		}
	}
	typedef map<string, vector<segment *> >::iterator it_type;
	for (it_type a = A.begin(); a!=A.end(); a++){
		for (int i =0; i < a->second.size(); i++){
			segments.push_back(a->second[i]);
		}
	}
	return segments;
}




//================================================================================================
//write out to file functions

void write_out(string FILE, map<string, interval_tree *> A){
	typedef map<string, interval_tree *>::iterator it_type;
	int N;
	ofstream FHW;
	
	FHW.open(FILE);
	
	for(it_type c = A.begin(); c != A.end(); c++) {
		N 	= A[c->first]->getDepth();
		merged_interval * Array  = new merged_interval[N];
		A[c->first]->insert_into_array(Array,N);
		for (int i = 0; i < N;i++){
			for (int j = 0; j < Array[i].intervals.size(); j++ ){
				if (Array[i].intervals[j].forward_x.size()>0 and Array[i].intervals[j].reverse_x.size()>0){
					FHW<<"#"<<c->first<<","<<to_string(Array[i].intervals[j].start)<<","<<to_string(Array[i].intervals[j].stop)<<endl;
					FHW<<"~forward"<<endl;
					for (int k = 0; k < Array[i].intervals[j].forward_x.size(); k++ ){
						FHW<<to_string(int(Array[i].intervals[j].forward_x[k]))<<","<<to_string(Array[i].intervals[j].forward_y[k])<<endl;
					}
					FHW<<"~reverse"<<endl;
					for (int k = 0; k < Array[i].intervals[j].reverse_x.size(); k++ ){
						FHW<<to_string(int(Array[i].intervals[j].reverse_x[k]))<<","<<to_string(Array[i].intervals[j].reverse_y[k])<<endl;
					}
				}
			}
		}
		delete[] Array;
	}
}

void write_out_bidir_fits( vector<segment*> segments, 
	map<int, map<int, bidir_preds> > G, params * P){

	int N 			= segments.size();
	string out_dir 	= P->p["-o"]; 
	string out_file_template 	= out_dir+"model_fits_out_bidir_" + P->p["-chr"]+"_";
	string out_file 			= check_file(out_file_template, 1);
	
	ofstream FHW;
	
	FHW.open(out_file_template+ "1" + ".bed");
	//FHW<<get_header(P);
	double scale 	= stod(P->p["-ns"]);

	typedef map<int, map<int, bidir_preds> >::iterator it_type_2;
	typedef map<int, bidir_preds>::iterator it_type_3;
	double center, left , right; 
	for (it_type_2 s = G.begin(); s!=G.end(); s++){
		for (it_type_3 b = s->second.begin(); b!=s->second.end(); b++ ){
			if (b->second.argK!=0){
				for (int i = 0; i < b->second.arg_all_model_components.all_components.size(); i++){
					center 	= b->second.arg_all_model_components.all_components[i].mu*scale + segments[s->first]->start;
					left 	= center-b->second.arg_all_model_components.all_components[i].si*scale;
					right 	= center+b->second.arg_all_model_components.all_components[i].si*scale;	
					FHW<<( 	segments[s->first]->chrom + "\t" + 
							to_string(int(left)) + "\t" + 
							to_string(int(right))+" \n");
					cout<<( 	segments[s->first]->chrom + "\t" + 
							to_string(int(left)) + "\t" + 
							to_string(int(right))+" \n")	;
				}
			}

		}
	}
}

void write_out_bidirs(map<string , vector<vector<double> > > G, string out_dir, 
	string job_name,int job_ID, params * P){
	typedef map<string , vector<vector<double> > >::iterator it_type;
	ofstream FHW;
	FHW.open(out_dir+ job_name+ "-" + to_string(job_ID)+ "_prelim_bidir_hits.bed");
	FHW<<P->get_header(4);
	int ID 	= 0;
	for (it_type c = G.begin(); c!=G.end(); c++){
		for (int i = 0; i < c->second.size(); i++){
			FHW<<c->first<<"\t"<<to_string(int(c->second[i][0]))<<"\t"<<to_string(int(c->second[i][1]))<<"\tME_"<<to_string(ID)<<endl;
			ID++;
		}
	}
	FHW.close();
}

void append_noise_intervals(string noise_bed_file, string bidir_file){
	ofstream bidir_pred;	
	ifstream FH(noise_bed_file);
  	random_device rd;
	mt19937 mt(rd());
	
	uniform_real_distribution<double> dist(0, 1);
	
  	bidir_pred.open (bidir_file,  ofstream::out |  ofstream::app);
	if (FH){
		string chrom, line;
		vector<string> lineArray;
		double start, stop,l, center;
		while(getline(FH, line)){
			lineArray=splitter(line, "\t");
			chrom 	= lineArray[0];
			start, stop 	= stod(lineArray[1]), stod(lineArray[2]);
			l 				= stop - start;
			center 	= (stop+start)/2.;
			if (chrom.size()==4  and l > 10000 and dist(mt)<0.05 ){
				bidir_pred<<chrom+"\t"+to_string(int(center-1000))	+"\t"+ to_string(int(center+1000))+"\tNOISE"<<endl;
			}
		}
	}else{
		printf("COULDNT OPEN NOISE FILE\n");
	}
	FH.close();
	bidir_pred.close();
}






string getp4_param_header(params * P){
	string header 	= "";
	return header;
}

void write_out_MLE_model_info(vector<final_model_output> A, params * P, string job_name, int JOB_ID ){
	string out_file_dir 	= P->p4["-o"];
	//so we want to write out two files:
	//(1) a specific file format with all parameter estimates listed
	//(2) a bed file showing final bidir predictions (1 standard deviation)

	ofstream FHW_bed;
	ofstream FHW_config;
	
	FHW_bed.open(out_file_dir+ job_name+ "-" +to_string(JOB_ID)+ "_bidirectional_hits_intervals.bed");
	FHW_bed<<P->get_header(4);
	for (int i = 0; i < A.size(); i++){
		FHW_bed<<A[i].write_out_bed();
	}

	FHW_bed.close();
	FHW_config.close();
}





vector<segment*> load_intervals_of_interest(string FILE, map<int, string>&  IDS, int pad,string spec_chrom){
	ifstream FH(FILE);
	vector<segment *> G;
	int ct 	= 1;
	if (FH){
		string line, chrom;
		int start, stop;
		int 	i = 0;
		vector<string>lineArray;
		string strand; 
		while(getline(FH, line)){
			lineArray=splitter(line, "\t");
			if (lineArray[0].substr(0,1)!="#" and lineArray.size()==4){
				if (lineArray.size() > 3){
					IDS[i] 		= lineArray[3];
				}
				if (lineArray.size() > 4){
					strand 		= lineArray[4];
				}else{
					strand 		= ".";
				}
				chrom=lineArray[0], start=max(stoi(lineArray[1])-pad, 0), stop=stoi(lineArray[2]) + pad;
				if (spec_chrom=="all" or spec_chrom==chrom){
					segment * S 	= new segment(chrom, start, stop,i,strand);
					G.push_back(S);
				}
				i++;
			}
		}
	}else{
		printf("couldn't open %s for reading\n", FILE.c_str() );
	}
	return G;
}




vector<segment *>  combind_bidir_fits_with_intervals_of_interest(vector<final_model_output> A, vector<segment *> FSI ){
	map<string, vector< vector<double> >> a;
	for (int i = 0; i < A.size(); i++){
		vector<vector<double>> d 	= A[i].get_bounds();
		for (int u = 0 ; u < d.size(); u++ ){
			d[u].push_back(0);
			a[A[i].chrom].push_back(d[u]);
		}
		
	}
	map<string, vector<segment *> > fsi;
	for (int i = 0; i < FSI.size(); i++){
		fsi[FSI[i]->chrom].push_back(FSI[i]);
	}
	vector<segment *> final_segments;
	typedef map<string, vector<segment *> >::iterator it_type;
	typedef vector<segment *>::iterator it_type_2;
	typedef map<string, vector< vector<double> >> ::iterator it_type_3;
	
	//want to check if they are sorted?
	for (it_type i = fsi.begin(); i !=fsi.end(); i++){
		for (int j = 1; j < i->second.size(); j++){
			if (i->second[j-1]->start > i->second[j]->start){
				printf("ERROR:, elongation support bed file is not sorted...\n");
			}
		}
	}
	for (it_type_3 i = a.begin(); i !=a.end(); i++){
		a[i->first] = interval_sort(i->second);
	}
	int N,j;
	for (it_type i = fsi.begin(); i !=fsi.end(); i++){
		if (a.find(i->first) != a.end()){
			j=0,N=a[i->first].size();
			vector<vector<double> > current 	= a[i->first];
			for (it_type_2 s 	= i->second.begin(); s != i->second.end(); s++ ){
				while (j < N and current[j][1] < (*s)->start){
					j++;
				}
				while (j < N and current[j][0] < (*s)->stop ){
					(*s)->add_fitted_bidir(current[j]);
					current[j][current[j].size()-1] 	= 1;
					j++;
				}
			}
			int p=0;
			for (int u = 0; u < current.size(); u++){
				if (current[u][current[u].size()-1] == 0){
					p++;
					segment * NS 	= new segment(i->first, current[u][0]-3000, current[u][1]+3000, -p );
					NS->add_fitted_bidir(current[u]);
					i->second.push_back(NS);
				}
			}
		}
	}	
	for (it_type i = fsi.begin(); i !=fsi.end(); i++){
		
		for (int j = 0; j < i->second.size(); j++){
			if (!i->second[j]->fitted_bidirs.empty()){
				final_segments.push_back(i->second[j]);
			}

		}
	}
	return final_segments;


}

void write_config_file_model_fits(vector<final_model_output> A, map<int, string> IDS, params * P ){
	double scale 	= stod(P->p4["-ns"]);
	string out_dir 	= P->p4["-o"];
	ofstream FHW_config;
	FHW_config.open(out_dir+"model_fits.txt");
	typedef vector<final_model_output>::iterator it_type;
	typedef vector<rsimple_c>::iterator it_type_2;
	string chrom_ID, ID;
	string  sigmas, mus, lambdas, pis_N, weights_N, weights_F, weights_R, bounds_F, bounds_R, pis_F, pis_R;
	string ll,N;
	int k 	= 0;
	FHW_config<<P->get_header(4);
	FHW_config<<"#chrom:start-stop\tunique identifier\tlog-likelihood,N\tLoading Positions\tVariance in Loading\tInitiation Lengths\tStrand Prob. Loading\tWeights Loading\tWeights Forward\tWeights Reverse\tForward Support\tReverse Support\tForward Strand Prob\tReverse Strand Prob\n";
	for (it_type a = A.begin(); a!= A.end(); a++){
		sigmas= "", mus= "", lambdas="", pis_N= "", weights_N= "", weights_F= "", weights_R= "", bounds_F= "", bounds_R= "", pis_F="", pis_R="";
		chrom_ID="", ID="", ll="";
		chrom_ID=(*a).chrom + ":"+ to_string((*a).start) + "-" + to_string((*a).stop);
		ID 			= ".";
		if (IDS.find((*a).ID)  != IDS.end() ){
			ID 		= IDS[(*a).ID];
		}
		string comma 	= "";
		k 	= 0;
		ll 	= to_string((*a).components[0].ps[1]);
		N 	= to_string((*a).components[0].ps[13]);
		for (it_type_2 r =(*a).components.begin(); r!= (*a).components.end(); r++ ){
			if (k < (*a).components.size() -1 ){
				comma = ",";
			}else{
				comma = "";
			}
			mus+=to_string((*r).ps[2]*scale + (*a).start) + comma;
			sigmas+=to_string((*r).ps[3]*scale ) + comma;
			lambdas+=to_string((1. / (*r).ps[4])*scale) + comma;
			pis_N+=to_string( (*r).ps[6] ) + comma;
			
			weights_N+=to_string((*r).ps[5] ) + comma;
			weights_F+=to_string((*r).ps[7] ) + comma;
			weights_R+=to_string((*r).ps[8] ) + comma;
			bounds_F+=to_string((*r).ps[9]*scale + (*a).start ) + comma;
			bounds_R+=to_string((*r).ps[10]*scale + (*a).start )+ comma;
			pis_F+=to_string((*r).ps[11]) + comma;
			pis_R+=to_string((*r).ps[12]) + comma;
			k++;
		}
		FHW_config<<chrom_ID+"\t"+ID+"\t"+ ll + ","+ N + "\t" + mus+"\t" + sigmas+"\t" + lambdas+ "\t" + pis_N+ "\t";
		FHW_config<<weights_N+ "\t" + weights_F+ "\t" + weights_R+ "\t" + bounds_F+ "\t" + bounds_R + "\t";
		FHW_config<<pis_F+ "\t" + pis_R+ "\n" ;
	}



}


void write_gtf_file_model_fits(vector<final_model_output> A, params * P){
	string out_dir 	= P->p4["-o"];

	ofstream FHW_gff;
	FHW_gff.open(out_dir+"model_fits.gff3");
	
	typedef vector<final_model_output>::iterator it_type;
	typedef vector<rsimple_c>::iterator it_type_2;
	string chrom;
	int center,std;
	vector<rsimple_c> components;
	double scale 	= stod(P->p4["-ns"]);
	int start, stop, left, right;
	string line;
	string parent,bidir,L,R,MR; 
	int 	i 		= 0;
	for (it_type a = A.begin(); a!= A.end(); a++){
		chrom 		= (*a).chrom;
		components 	= (*a).components;
		start 		= (*a).start;
		for (it_type_2 r =(*a).components.begin(); r!= (*a).components.end(); r++ ){
			//human15.1 . gene            214301  215772 . +   . ID=HsG8283
			left 	= (*r).ps[10]*scale + start;
			right 	= (*r).ps[9]*scale + start;
			
			center 	= start + (*r).ps[2]*scale;
			std 	= (*r).ps[3]*scale*0.5 + (scale/(*r).ps[4] ) ;
			parent 	= chrom +"\t.\tgene\t";
			parent+= to_string(center) + "\t" + to_string(right) +"\t.\t.\t.\tID=Bidir_"+to_string(i);
			
			
			bidir 	= chrom +"\tEMGU\tCDS\t";
			bidir 	+= to_string(center) + "\t" + to_string(center+std) + "\t.\t+\t.\tParent=Bidir_"+to_string(i);
			bidir 	+=";Name=Center "+ to_string(center) +  ", std " + to_string(int((*r).ps[3]*scale  )) + ", init " + to_string(int((scale/(*r).ps[4] ) ));
			FHW_gff<<parent<<endl;
			FHW_gff<<bidir<<endl;
			if (((*r).ps[7] > 0.001)){
				R 		= chrom +"\tEMGU\tCDS\t";
				R 		+= to_string(right) + "\t" + to_string(std+right) + "\t.\t+\t.\tParent=Bidir_"+to_string(i);
				R 		+=";Name=Forward Strand Elongation Component";
				FHW_gff<<R<<endl;
			}
			
			if ((*r).ps[8] > 0.001){
				L 		= chrom +"\tEMGU\tCDS\t";
				L 		+= to_string(left-std) + "\t" + to_string(left) + "\t.\t-\t.\tParent=Bidir_"+to_string(i);
				L+=";Name=Reverse Strand Elongation Component";
				FHW_gff<<L<<endl;
			}
			
			//FHW_gff<<L<<endl;
			
			line="";
			i++;
		}

	}
	FHW_gff.close();
}


vector<segment* > insert_bedgraph_to_segment_single(map<string, vector<segment *> > A , string FILE, int rank){
	
	//want instead to make this interval tree
	ofstream FHW;
	map<string, vector<merged_interval*> > merged_FSI 	= segments_to_merged_intervals( A, FHW);
	map<string, interval_tree *> AT;
	typedef map<string, vector<merged_interval*>>::iterator it_type_4;
	for(it_type_4 c = merged_FSI.begin(); c != merged_FSI.end(); c++) {
		AT[c->first] 	= new interval_tree();
		AT[c->first]->construct(c->second);
	}
	

	int start, stop, N, j;
	double coverage;
	N 	= 0,j 	= 0;
	int strand 	= 1;
	int o_st, o_sp;
	vector<string> lineArray;
	string chrom, prevchrom, line;
	vector<segment *> segments;
	double center;
	ifstream FH(FILE );
	if (FH){
		prevchrom="";
		while (getline(FH, line)){
			lineArray 	= splitter(line, "\t");
			chrom 		= lineArray[0];
			start=stoi(lineArray[1]),stop=stoi(lineArray[2]), coverage = stod(lineArray[3]);
			center 	= (stop + start) /2.;
			if (AT.find(chrom)!=AT.end()){
				AT[chrom]->insert(center, coverage, strand);
			}			
		}
		FH.close();
		
	}else{
		cout<<"could not open: "<<FILE<<endl;
		segments.clear();
		return segments;
	}

	//now we want to get all the intervals and make a vector<segment *> again...
	typedef 	map<string, vector<segment *> >::iterator it_type_5;
	for(it_type_5 c = A.begin(); c != A.end(); c++) {
		for (int i = 0; i < c->second.size(); i++){
			center 	= (c->second[i]->stop + c->second[i]->start) /2.;
			interval FOUND 	= AT[c->first]->get_interval(c->second[i]->start, c->second[i]->stop );
			if (not FOUND.EMPTY and FOUND.forward_x.size()){
				segment * S 	= new segment(c->first, FOUND.start, FOUND.stop, FOUND.ID);
				S->counts 		= FOUND.counts;
				S->parameters 	= FOUND.parameters;
				for (int u = 0 ; u < FOUND.forward_x.size(); u++){
					vector<double> curr(2);
					curr[0] 	= FOUND.forward_x[u], curr[1] 	= FOUND.forward_y[u];
					
					S->forward.push_back(curr);
				}
				segments.push_back(S);

			}

		}
	}
	return segments;
}

vector<segment* > insert_bedgraph_to_segment_joint(map<string, vector<segment *> > A , 
	string forward, string reverse, int rank, ofstream& FHW){
	
	//want instead to make this interval tree
	FHW<<"(load) making merged intervals: "<<A.size()<<endl;
	FHW.flush();
	
	map<string, vector<merged_interval*> > merged_FSI 	= segments_to_merged_intervals( A, FHW);
	FHW<<"(load) made merged intervals"<<endl;
	FHW.flush();
	map<string, interval_tree *> AT;
	typedef map<string, vector<merged_interval*>>::iterator it_type_4;
	for(it_type_4 c = merged_FSI.begin(); c != merged_FSI.end(); c++) {
		AT[c->first] 	= new interval_tree();
		AT[c->first]->construct(c->second);
	}
	FHW<<"(load) interval tree"<<endl;
	FHW.flush();
	

	int start, stop, N, j;
	double coverage;
	N 	= 0,j 	= 0;
	int strand;
	int o_st, o_sp;
	vector<string> lineArray;
	string chrom, prevchrom, line;
	vector<segment *> segments;
	double center;
	vector<string> FILES(2);
	FILES[0]=forward, FILES[1]=reverse;
	string FILE;
	for (int i =0; i < 2; i++){
		if (i==0){
			strand=1;
		}else{
			strand=-1;
		}
		FILE=FILES[i];
		ifstream FH(FILE);

		if (FH){
			prevchrom="";
			while (getline(FH, line)){
				lineArray 	= splitter(line, "\t");
				chrom 		= lineArray[0];
				start=stoi(lineArray[1]),stop=stoi(lineArray[2]), coverage = abs(stod(lineArray[3]));
				center 	= (stop + start) /2.;
				if (AT.find(chrom)!=AT.end()){
					AT[chrom]->insert(center, coverage, strand);
				}			
			}
			FH.close();
			
		}else{
			cout<<"could not open forward bedgraph file: "<<FILE<<endl;
			segments.clear();
			return segments;
		}
	}

	FHW<<"read in file"<<endl;
	FHW.flush();
	//now we want to get all the intervals and make a vector<segment *> again...
	typedef 	map<string, vector<segment *> >::iterator it_type_5;
	for(it_type_5 c = A.begin(); c != A.end(); c++) {
		for (int i = 0; i < c->second.size(); i++){
			center 	= (c->second[i]->stop + c->second[i]->start) /2.;
			interval FOUND 	= AT[c->first]->get_interval(int(c->second[i]->start), int(c->second[i]->stop) );
			if (not FOUND.EMPTY and FOUND.forward_x.size()){
				segment * S 	= new segment(c->first, FOUND.start, FOUND.stop, FOUND.ID, FOUND.STRAND);
				S->parameters 	= FOUND.parameters;
				for (int u = 0 ; u < FOUND.forward_x.size(); u++){
					vector<double> curr(2);
					curr[0] 	= FOUND.forward_x[u], curr[1] 	= FOUND.forward_y[u];
					
					S->forward.push_back(curr);
				}
				for (int u = 0; u < FOUND.reverse_x.size(); u++){
					vector<double> curr(2);
					curr[0] 	= FOUND.reverse_x[u], curr[1] 	= FOUND.reverse_y[u];
					S->reverse.push_back(curr);
					
				}
				segments.push_back(S);
			}
		}
	}
	FHW<<"made segments"<<endl;
	FHW.flush();
	return segments;
}

void write_out_single_simple_c(vector<single_simple_c> fits, map<int, string> IDS , params * P ){
	string out_dir 	= P->p5["-o"];
	double scale 	= stod(P->p5["-ns"]);
	ofstream FHW_gff;
	ofstream FHW_config;
	FHW_gff.open(out_dir+"single_model_fits.gff3");

	FHW_config.open(out_dir+"single_model_fits.txt");
	FHW_config<<P->get_header(5);
	FHW_config<<"#chrom:start-stop\tID/Name\tmu(loading position)\tsigma(variance in loading)\tmixture weights\tlog-likelihood\tData Point Number\n";
	
	string chrom, ID, strand, parent, loading, L, R;
	int start, stop, left, right;
	int center, std;
	string config_line;
	for (int i = 0; i < fits.size(); i++){
		chrom 	= fits[i].chrom;
		start 	= fits[i].st_sp[0];
		center 	= fits[i].ps[0]*scale  + start;
		std 		= fits[i].ps[1]*scale;
		ID 		= IDS[fits[i].st_sp[2]];
		if (ID.empty()){
			ID 	= to_string(i);
		}
		left 			= fits[i].ps[3]*scale + fits[i].st_sp[0], right 	= fits[i].ps[4]*scale + fits[i].st_sp[0];
		if (fits[i].ps[6] > pow(10,-2) and fits[i].ps[5] > pow(10,-2)){

			strand 	= ".";
		}else if (fits[i].ps[6] > pow(10,-2)){
			strand 	= "+";
		}else{
			strand 	= "-";
		}
		

		parent 	= chrom +"\tEMGU\tgene\t";
		parent+= to_string(fits[i].st_sp[0]) + "\t" + to_string(fits[i].st_sp[1]) +"\t.\t"+ strand+"\t.\tID="+ID + "\n";
		loading 	= "";
		L = "", R 	= "";		
			
		if (fits[i].ps[2] > pow(10,-1)){

			loading 	= chrom +"\tEMGU\tCDS\t";
			loading 	+= to_string(center-std) + "\t" + to_string(center+std) + "\t.\t" + strand+"\t.\tParent="+ID;
			loading 	+=";Name=Center "+ to_string(center) +  ", std " + to_string( std) + "\n";
			if (strand == "."){
				L 		= chrom +"\tEMGU\tCDS\t";
				L 		+= to_string(left-std) + "\t" + to_string(left) + "\t.\t.\t.\tParent="+ID;
				L+=";Name=Reverse Strand Elongation Component\n";

				R 		= chrom +"\tEMGU\tCDS\t";
				R 		+= to_string(right) + "\t" + to_string(std+right) + "\t.\t.\t.\tParent="+ID;
				R 		+=";Name=Forward Strand Elongation Component\n";	
			}else if (strand == "+"){
				R 		= chrom +"\tEMGU\tCDS\t";
				R 		+= to_string(right) + "\t" + to_string(std+right) + "\t.\t+\t.\tParent="+ID;
				R 		+=";Name=Forward Strand Elongation Component\n";	
			}else{
				L 		= chrom +"\tEMGU\tCDS\t";
				L 		+= to_string(left-std) + "\t" + to_string(left) + "\t.\t-\t.\tParent="+ID;
				L 		+=";Name=Forward Strand Elongation Component\n";	
			}
		}else{
			strand 	= ".";
			parent 	= chrom +"\tEMGU\tgene\t";
			parent+= to_string(fits[i].st_sp[0]) + "\t" + to_string(fits[i].st_sp[1]) +"\t.\t"+ strand+"\t.\tID="+ID + "\n";
			L 		= chrom +"\tEMGU\tCDS\t";
			L 		+= to_string(left-std) + "\t" + to_string(left) + "\t.\t.\t.\tParent="+ID;
			L 		+=";Name=Elongation Only\n";	
			R 		= chrom +"\tEMGU\tCDS\t";
			R 		+= to_string(right) + "\t" + to_string(std+right) + "\t.\t.\t.\tParent="+ID;
			R 		+=";Name=Elongation Only\n";	
			loading 	= "";
		
		}
			
		
		FHW_gff<<parent;
		FHW_gff<<loading;
		FHW_gff<<R;
		FHW_gff<<L;

		//write_out_config_file
		config_line 	= chrom +":" + to_string(fits[i].st_sp[0]) + "-" + to_string(fits[i].st_sp[1]) + "\t";
		config_line    += ID + "\t";
		config_line    += to_string((fits[i].ps[0]*scale  + start)) + "\t" + to_string(fits[i].ps[1]*scale) + "\t" + to_string(fits[i].ps[2]) + "," + to_string(fits[i].ps[5]) + "," + to_string(fits[i].ps[6]);
		config_line    += "\t" + to_string(fits[i].ps[7]) +  "\t" +to_string(fits[i].ps[8]) + "\n";
		FHW_config<<config_line;
		



	}

}

void write_out_single_simple_c_marks(vector<single_simple_c> fits, map<int, string> IDS , params * P ){
	string out_dir 	= P->p5["-o"];
	double scale 	= stod(P->p5["-ns"]);
	ofstream FHW;
	
	FHW.open(out_dir+"refined_peaks.bed");
	FHW<<P->get_header(5);
	string chrom;
	int start, stop;
	double center, sig, w;
	string line, params;
	for (int i = 0; i < fits.size(); i++){
		if (fits[i].ps[0]*scale > 0 and fits[i].ps[8] > nINF  ){
			chrom 	= fits[i].chrom;
			center 	= fits[i].ps[0]*scale  + fits[i].st_sp[0];
			sig 	= fits[i].ps[1]*scale;
			start 	= center - sig;
			stop 	= center + sig;

			params 	= to_string(fits[i].ps[0]*scale) + "," + to_string(fits[i].ps[1]*scale) + "," + to_string(scale/ fits[i].ps[2] );
			for (int u = 3; u < 10; u++){
			
				if (u == 5){
					params+=(","+to_string(fits[i].ps[u]*scale));
				}else{
					params+=(","+to_string(fits[i].ps[u]));
				}
			}
			line 	= chrom + "\t" + to_string(start) + "\t" + to_string(stop) + "\t" + params + "\n";
			FHW<<line;
		}
	}
	FHW.close();
}


void collect_all_tmp_files(string dir, string job_name, int nprocs, int job_ID){
	int c 	= 0;
	time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
    // for more information about date/time format
    strftime(buf, sizeof(buf), "%m_%d_%H_%M", &tstruct);
    string DT 	= buf;
	string OUT 		= dir+ job_name + "-" + to_string(job_ID) +"_" + DT+ ".log";
	ofstream FHW(OUT);
	for (int rank = 0; rank < nprocs; rank++){
		string FILE 	= dir+"tmp_" + job_name+ "-" +to_string(job_ID) + "_" + to_string(rank) + ".log";
		string line;
		ifstream FH(FILE);
		if (FH){
			if (rank!=0){
				FHW<<"=======MPI Call: " + to_string(rank) +"=======\n";
			}
			while (getline(FH, line)){
				if ("#" != line.substr(0,1) or rank==0){
					FHW<<line<<endl;
				}
			}

			FH.close();
			remove( FILE.c_str()) ;
			c++;
		}else{
			printf("HERE?????");
		}
	}
}


void get_noise_mean_var(string noise_file, string bedgraph, double * mean, double * var){

	map<string, vector<vector<double> >> G;
	ifstream FH(noise_file);
	if (FH){
		string line, chrom, start, stop;
		vector<string> lineArray;

		while (getline(FH, line)){
			lineArray 	= splitter(line, "\t");
			chrom=lineArray[0], start=lineArray[1], stop=lineArray[2];
			vector<double> current(3);
			current[0]=stof(start), current[1]=stof(stop), current[2]=0;

			G[chrom].push_back(current);
		
		
		}
	}else{
		cout<<"Coudn't open: "<<noise_file<<endl;
	}

	FH.close();
	ifstream BED(bedgraph);
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
			while (j < N and G[chrom][j][1]<start){
				j++;
			}
			if (j < N and G[chrom][j][0]<stop){
				G[chrom][j][2]+=(stop-start)*cov;
			}
			prevchrom=chrom;
		}
		
	}else{
		cout<<"Coudn't open: "<<bedgraph<<endl;

	}
	BED.close();
	typedef map<string, vector<vector<double>> >::iterator it_type;
	double SUM=0, N=0;
	for (it_type g=G.begin(); g!=G.end(); g++){
		for (int i = 0; i < g->second.size(); i++){
			if (g->second[i][2] > 0){
				SUM+=(g->second[i][2]/(g->second[i][1]-g->second[i][0]));
				N+=1;
			}

		}
	}
	*(mean) 	= SUM/N;
	double VAR=0;
	for (it_type g=G.begin(); g!=G.end(); g++){
		for (int i = 0; i < g->second.size(); i++){
			if (g->second[i][2] > 0){
				VAR+=pow((g->second[i][2]/(g->second[i][1]-g->second[i][0] ))   -(*mean), 2);
			}

		}
	}
	
	*(var) 	= VAR/N;


}


void write_out_models_from_free_mode(map<int, map<int, vector<simple_c_free_mode>  > > G, 
	params * P, int job_ID,map<int, string> IDS){
	//========================================================================================
	//write out each model parameter estimates
	double scale 	= stof(P->p["-ns"]);
	double penality = stof(P->p["-ms_pen"]) ;
	string out_dir 	= P->p["-o"];
	ofstream FHW;


	FHW.open(out_dir+  P->p["-N"] + "-" + to_string(job_ID)+  "_K_models_MLE.tsv");
	FHW<<P->get_header(0);
	
	typedef map<int, map<int, vector<simple_c_free_mode>  > >::iterator it_type_1;
	typedef map<int, vector<simple_c_free_mode>  > ::iterator it_type_2;
	typedef vector<simple_c_free_mode>::iterator it_type_3;
	
	typedef map<int, string>::iterator it_type_IDS;
	int IN=0;
	string mus="", sis="", ls="", wEMs="", wPIs="",forward_bs="", forward_ws="",forward_PIs="",reverse_as="", reverse_ws="",reverse_PIs="";
	double w_thresh 	= 0.;
	double ALPHA_2 	= stof(P->p["-ALPHA_2"]);
	string mu, sigma, lambda, pos, neg, ll, pi, w,ra,fb,rw,fw,fp;
	int start;
	string chrom;

	string INFO 	= "";
			
	FHW<<"#ID|chromosome:start-stop|forward strand coverage, reverse strand coverage"<<endl;
	FHW<<"#model complexity,log-likelihood"<<endl;
	FHW<<"#mu_k\tsigma_k\tlambda_k\tpi_k\tfp_k\tw_[p,k],w_[f,k],w_[r,k]\tb_[f,k]\ta_[r,k]"<<endl;
	
	for (it_type_1 s = G.begin(); s!=G.end(); s++){ //iterate over each segment
		FHW<<">" + IDS[s->first]+ "|";
		for (it_type_2 k 	= s->second.begin(); k != s->second.end(); k++){//iterate over each model_complexity
			for (it_type_3 c = k->second.begin(); c!=k->second.end(); c++){
				chrom 		= (*c).chrom;
				INFO 		= chrom + ":" + to_string((*c).ID[1])+"-"+to_string((*c).ID[2]);
				pos 		= to_string((*c).SS[1]);
				neg 		= to_string((*c).SS[2]);
				
			}
		}
		FHW<<INFO<<"|"<<pos+","+neg<<endl;


		for (it_type_2 k 	= s->second.begin(); k != s->second.end(); k++){//iterate over each model_complexity
			string mus 		= "", sigmas="", lambdas="",pis="", ws= ""  ,fbs="",ras="", fws="", rws="", fps="";
			string pos 		= "", neg="";
			string k_header = "~"+ to_string(k->first)+",";
			int NN 			= k->second.size();
			int ii 			= 0;
			for (it_type_3 c = k->second.begin(); c!=k->second.end(); c++){
				chrom 		= (*c).chrom;
				start 		= (*c).ID[1];
				mu 			= to_string((*c).ps[0]*scale + (*c).ID[1] );
				sigma 		= to_string((*c).ps[1]*scale);
				lambda 		= to_string(scale/(*c).ps[2]);
				pi 			= to_string( (*c).ps[4]);
				w  			= to_string( (*c).ps[3]);
				fw 			= to_string( (*c).ps[6]); 
				rw 			= to_string( (*c).ps[9]);
				ra 			= to_string( scale*(*c).ps[8]   + (*c).ID[1] );
				fb 			= to_string( scale*(*c).ps[5]  + (*c).ID[1]);
				fp 			= to_string( scale*(*c).ps[11]  );
				ll 			= to_string((*c).SS[0]);
				if (ii +1 < NN){


					mus+=mu+",";
					sigmas+=sigma+",";
					lambdas+=lambda+",";
					pis+=pi+",";
					ws+=w+ "," + fw+ "," + rw + "|";
					ras+=ra+",";
					fbs+=fb+",";
					fps+=fp+",";
				}else{
					mus+=mu ;
					sigmas+=sigma ;
					lambdas+=lambda ;
					pis+=pi ;
					ws+=w+ "," + fw+ "," + rw ;	
					fbs+=fb;
					ras+=ra;
					fps+=fp ;
				}
				ii++;
			}
			k_header 		+=ll+ "\t";
			FHW<<k_header;
		
			if (k->first>0){
				FHW<<mus+"\t"+sigmas+"\t"+lambdas+"\t" + pis+"\t" + fps+ "\t" + ws + "\t" + fbs+"\t" +ras ;
			}
			FHW<<endl;
		}
	}

	//========================================================================================
	//perform model selection, write out bed files
	ofstream FHW_bed;


	FHW_bed.open(out_dir+  P->p["-N"] + "-" + to_string(job_ID)+  "_bidirectional_hits_intervals.bed");
	FHW_bed<<P->get_header(0);
	double Penality 	= 100;
	double BIC_score;
	string M,S,L,PI, W,FW, RW, FP, LL, POS,NEG;
	for (it_type_1 s = G.begin(); s!=G.end(); s++){ //iterate over each segment
		double BIC_MAX 		= INF;
		int arg_K 		= 0;
		for (it_type_2 k 	= s->second.begin(); k != s->second.end(); k++){//iterate over each model_complexity
			for (it_type_3 c = k->second.begin(); c!=k->second.end(); c++){
			
				BIC_score 	= -2*(*c).SS[0] + (k->first+1)*(Penality)*log((*c).SS[1] + (*c).SS[2]);
				if (BIC_score < BIC_MAX){
					arg_K 	= k->first;
					BIC_MAX = BIC_score;
				}
			}
		}
		if (arg_K!=0){
			for (it_type_3 c 	= s->second[arg_K].begin(); c!= s->second[arg_K].end(); c++   ){
				string chrom 		= (*c).chrom;
				int start 		= (*c).ID[1];
				double mu 			= ((*c).ps[0]*scale + (*c).ID[1] );
				double sigma 		= ((*c).ps[1]*scale);
				double lambda 		= scale/(*c).ps[2];
				if (sigma < 1000 and lambda < 1000 and (*c).ps[3] > 0.1 ){
					M  			= to_string((*c).ps[0]*scale + (*c).ID[1] );
					S 			= to_string((*c).ps[1]*scale);
					L 			= to_string(scale/(*c).ps[2]);
					PI 			= to_string( (*c).ps[4]);
					W 			= to_string( (*c).ps[3]);
					FW 			= to_string( (*c).ps[6]); 
					RW 			= to_string( (*c).ps[9]);
					FP 			= to_string( scale*(*c).ps[11]  );
					LL 			= to_string((*c).SS[0]);
					POS 		= to_string((*c).SS[1]);
					NEG 		= to_string((*c).SS[2]);

					INFO 		= IDS[s->first]+ "|" + M + "," + S + ","+ L +","+PI+","+W+","+FW+","+RW + ","+LL+ "," +POS+ ","+NEG  ;

					FHW_bed<<chrom+"\t"+to_string(int(mu-sigma-lambda))+"\t"+to_string(int(mu+sigma+lambda))<<"\t"<<INFO<<endl;
				}
			}
		}
	}
	




	

}

struct bidir_segment{
public:
	string chrom; 
	int start, stop, chrom_ID;
	string INFO;
	vector<string > INFOS;
	bidir_segment(){};
	bidir_segment(string c, int st, int sp, string info){
		chrom=c, start=st, stop=sp, INFO=info;
	}
	vector<vector<double> > get_parameters(){
		vector<string> lineArray;
		vector<vector<double>> centers; 
		for (int i = 0 ; i < INFOS.size(); i++){
			lineArray 	= splitter(INFOS[i], "_");
			vector<double> current(8);
			for (int i = 0; i < 8;i++){
				current[i] 	= stod(lineArray[i]);
			}
			centers.push_back(current);
		}
		return centers;

	}

};
vector<bidir_segment> bubble_sort_bidir_segment(vector<bidir_segment> X){
	bool changed=true;
	while (changed){
		changed=false;
		for (int i = 0; i < X.size()-1; i++  )	{
			if (X[i].start > X[i+1].start){ //sort by starting position
				bidir_segment copy 			= X[i];
				X[i] 					= X[i+1];
				X[i+1] 					= copy;
				changed=true;
			}
		}
	}
	return X;
}

vector<bidir_segment> get_merged(params * P){
	string FILE 		= P->p6["-i"];
	string spec_chrom 	= P->p6["-chr"];
	int pad 			= stoi(P->p6["-pad"]);
	ifstream FH(FILE);
	map<string, vector< bidir_segment >> G;
	int ct 	= 0;
	int i 	= 0;
	if (FH){
		string line, chrom;
		int start, stop;
		int 	i = 0;
		string strand; 
		vector<string> lineArray;
		while(getline(FH, line)){
			if (line.substr(0,1)!="#"){
				lineArray=splitter(line, "\t");
				
				chrom=lineArray[0], start=max(stoi(lineArray[1])-pad, 0), stop=stoi(lineArray[2])+pad;
				vector<int> current(3);
				current[0]=start-pad, current[1]=stop+pad, current[2]=i;

				if (chrom==spec_chrom or spec_chrom=="all"){
					G[chrom].push_back(bidir_segment(chrom, start, stop, lineArray[3]) );
					ct++;
				}
			}
			i++;
		}
	}else{
		printf("couldn't open %s for reading\n", FILE.c_str() );
	}

	//sort by starting
	typedef map<string, vector< bidir_segment  > >::iterator it_type;
	for (it_type g = G.begin(); g!=G.end(); g++){
		G[g->first] 	= bubble_sort_bidir_segment(g->second);
	}
	vector<bidir_segment > Merged;
	for (it_type g = G.begin(); g!=G.end(); g++ ){
		int j = 0, N = g->second.size();
		vector<bidir_segment> d 	= g->second;
		bidir_segment current;
		if (j < N){
			current 	= d[j];
		}
		while ( j < N ){
			while (j < N and  d[j].stop > current.start  and  d[j].start < current.stop   ){
				current.start 	= min(current.start, d[j].start);
				current.stop 	= max(current.stop, d[j].stop);
				current.INFOS.push_back(d[j].INFO);
				j+=1;
			}
			Merged.push_back(current);
			if (j < N){
				current 	= d[j];
			}	
		}
	}
	return Merged;
}

map<string, vector<segment *> > load_bidir_predictions(params * P, 
	vector<int> st_sp, map<string, int>& chrom_to_ID, 
		map<int, string>& ID_to_chrom ){
	int pad 	= stoi(P->p6["-pad"]);
	vector<bidir_segment> Merged 	= get_merged(P);
	//make 
	int i 	= 0;
	for (int m = 0 ; m < Merged.size(); m++){
		if (chrom_to_ID.find(Merged[m].chrom) == chrom_to_ID.end() ){
			chrom_to_ID[Merged[m].chrom] 	= i;
			ID_to_chrom[i]  	= Merged[m].chrom;
			i++;
		}
	}
				
	map<string, vector<segment *> > G;
	for (int i = st_sp[0]; i < st_sp[1]; i++){
		segment * S 	= new segment(Merged[i].chrom, Merged[i].start, Merged[i].stop);
		S->parameters 	= Merged[i].get_parameters();
		S->chrom_ID 	= chrom_to_ID[S->chrom];
		G[S->chrom].push_back(S);
	}
	

	return G;

}

vector<vector<int> > get_line_start_stops(params * P, int nprocs){
	vector<bidir_segment> Merged 	= get_merged(P);
	int counts 	= int(Merged.size()) / nprocs;
	vector<vector<int> > start_stop;
	int start = 0, stop 	= 0;
	for (int j = 0 ; j < nprocs; j++){
		if (j+1 < nprocs){
			start 	= j*counts;
			stop 	= (j+1)*counts;
		}else{
			start 	= j*counts;
			stop 	= Merged.size();
		}
		vector<int> st_sp(2);
		st_sp[0] 	= start, st_sp[1]=stop;
		start_stop.push_back(st_sp);
	}

	return start_stop;


}




void write_bootstrap(vector<boostrap_struct> bs, map<int, string> ID_to_chrom, params * P,int JOB_ID){
	string out_file_dir 	= P->p6["-o"];
	string job_name 		= P->p6["-N"];

	ofstream FHW_bed;
	
	FHW_bed.open(out_file_dir+ job_name+ "-" +to_string(JOB_ID)+ "_bootstrapped_bidirectional_hits_intervals.bed");
	FHW_bed<<P->get_header(6);
	for (int b =0 ; b < bs.size(); b++){
		FHW_bed<<bs[b].print_out(ID_to_chrom, P);
	}
	
}


vector<segment *> merge_intervals_of_interest(vector<segment *> IOI){
	vector<segment * > FSI;
	map<string,vector<segment * > > G;
	typedef map<string,vector<segment * > >::iterator it_type;
	for (int i = 0; i < IOI.size(); i++){
		G[IOI[i]->chrom].push_back(IOI[i]);
	}
	for (it_type c = G.begin(); c!=G.end();c++){
		G[c->first] = bubble_sort_segment(c->second);
	}
	for (it_type c = G.begin(); c!= G.end(); c++){
		int N 	= c->second.size();
		int j 	= 0;
		int start, stop,ct;
		while (j < N){
			start=c->second[j]->start,stop=c->second[j]->stop;
			ct = 0;
			vector<vector<double>> centers;
			while (j < N and stop > c->second[j]->start and start < c->second[j]->stop ){
				start = min(start, c->second[j]->start), stop = max(stop, c->second[j]->stop);
				vector<double> center(1);
				center[0]=(c->second[j]->start + c->second[j]->stop)/2.;
				centers.push_back(center);
				ct++;
				j++;
			}
			segment * current 	= new segment(c->first,start, stop);
			current->counts 	= ct;
			current->parameters = centers;

			FSI.push_back(current);
		}
	}
	return FSI;
	

}









