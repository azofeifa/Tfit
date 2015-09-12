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
}

segment::segment(){
	N 		= 0;
	counts 	= 1;
	XN 		= 0;
	ID 		= 0;
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

	//===================
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
	else if (stop < current->start and left!=NULL){
		return left->get_interval(start, stop);
	}
	else if (start > current->stop and right!=NULL){
		return right->get_interval(start, stop);
	}
	return empty;

}

//================================================================================================
//loading from file functions...need to clean this up...

vector<segment*> load_bedgraphs_total(string forward_strand, 
	string reverse_strand, int BINS, double scale, string spec_chrom){
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
	
	while (getline(FH, line)){
		lineArray=splitter(line, "\t");
		chrom=lineArray[0], start=stoi(lineArray[1]), stop=stoi(lineArray[2]), coverage=stod(lineArray[3]);
		if (chrom != prevChrom and (chrom==spec_chrom or spec_chrom=="all")  )  {
			FOUND 		= true;
			G[chrom] 	= new segment(chrom, start, stop );
		}
		if (FOUND and chrom!= spec_chrom and spec_chrom!= "all"){
			break;
		}
		G[chrom]->add2(1, double((stop + start) / 2.), coverage);
		prevChrom=chrom;

	}
	prevChrom="";
	FOUND=false;
	while (getline(FH2, line)){
		lineArray=splitter(line, "\t");
		chrom=lineArray[0], start=stoi(lineArray[1]), stop=stoi(lineArray[2]), coverage=stod(lineArray[3]);
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

map<string, vector<merged_interval*> > segments_to_merged_intervals(map<string, vector<segment *> > FSI){
	map<string, vector<interval>> G;	
	typedef map<string, vector<segment *> >::iterator it_type_2;
	for (it_type_2 i = FSI.begin(); i != FSI.end(); i++){
		for (int j = 0; j < i->second.size(); j++){
			G[i->first].push_back(interval(i->second[j]->chrom, i->second[j]->start, i->second[j]->stop, i->second[j]->ID)  );
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
	string forward_file, string reverse_file, int pad, string spec_chrom){
	
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
			B[i->first].push_back(S);
		}
	}

	int N,j;
	//we want to merge the overlaping calls
	for (it_type_2 i = B.begin(); i!=B.end();i++){
		N 	= i->second.size(), j 	= 1;
		S 	= i->second[0];
		S->centers.push_back((S->start + S->stop) /2.);
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
				start=stoi(lineArray[1]),stop=stoi(lineArray[2]), coverage = stod(lineArray[3]);
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
				start=stoi(lineArray[1]),stop=stoi(lineArray[2]), coverage = stod(lineArray[3]);

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

void write_out_bidirs(map<string , vector<vector<double> > > G, string out_dir){
	typedef map<string , vector<vector<double> > >::iterator it_type;
	ofstream FHW;
	FHW.open(out_dir+"prelim_bidir_hits.bed");
	for (it_type c = G.begin(); c!=G.end(); c++){
		for (int i = 0; i < c->second.size(); i++){
			FHW<<c->first<<"\t"<<to_string(int(c->second[i][0]))<<"\t"<<to_string(int(c->second[i][1]))<<endl;
		}
	}
}

string getp4_param_header(params * P){
	string header 	= "";
	return header;
}

void write_out_MLE_model_info(vector<final_model_output> A, params * P ){
	string out_file_dir 	= P->p4["-o"];
	//so we want to write out two files:
	//(1) a specific file format with all parameter estimates listed
	//(2) a bed file showing final bidir predictions (1 standard deviation)

	ofstream FHW_bed;
	ofstream FHW_config;
	
	FHW_bed.open(out_file_dir+"bidirectional_hits_intervals.bed");
	FHW_config.open(out_file_dir+"parameters_out.tsv");
	
	for (int i = 0; i < A.size(); i++){
		FHW_bed<<A[i].write_out_bed();
		FHW_config<<A[i].write_out_config();
	}

	FHW_bed.close();
	FHW_config.close();
}





vector<segment*> load_intervals_of_interest(string FILE, map<int, string>&  IDS, int pad){
	ifstream FH(FILE);
	vector<segment *> G;
	int ct 	= 1;
	if (FH){
		string line, chrom;
		int start, stop;
		int 	i = 0;
		vector<string>lineArray;
		while(getline(FH, line)){
			lineArray=splitter(line, "\t");
			if (lineArray.size() > 3){
				IDS[i] 		= lineArray[3];
			}
			chrom=lineArray[0], start=max(stoi(lineArray[1])-pad, 0), stop=stoi(lineArray[2]);
			segment * S 	= new segment(chrom, start, stop,i);
			G.push_back(S);
			i++;
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
			for (int u = 0; u < current.size(); u++){
				if (current[u][current[u].size()-1] == 0){
					segment * NS 	= new segment(i->first, current[u][0]-3000, current[u][1]+3000 );
					NS->add_fitted_bidir(current[u]);
					i->second.push_back(NS);
				}
			}
		}
	}	
	for (it_type i = fsi.begin(); i !=fsi.end(); i++){
		for (int j = 0; j < i->second.size(); j++){
			final_segments.push_back(i->second[j]);

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
		ID 		= IDS[(*a).ID];
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
	map<string, vector<merged_interval*> > merged_FSI 	= segments_to_merged_intervals( A);
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











