#include <string>
#include <vector>
#include "load.h"
#include "split.h"
#include <iostream>
#include <fstream>
#include <map>
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
	//we also want to get rid of those data points that we don't need
	//i.e. the ones where there is no data coverage values on either the 
	//forward or reverse strands
	int realN 		= 0;
	for (int i = 0; i < BINS;i++){
		if (X[1][i]>0 or X[2][i]>0){
			realN++;
		}
	}

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
	//clear previous memory
	for (int i = 0; i < 3; i ++){
		delete X[i];
	}
	delete X;
	X 				= newX;
	XN 				= realN;
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
		if (segments[i]->forward.size() > 0 or segments[i]->reverse.size() > 0 ){
			segments[i]->bin(BINS, scale);
		}
	}
}

interval::interval(){};
interval::interval(string chr, int st, int sp){
	chrom 	= chr;
	start 	= st;
	stop 	= sp;
};

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


//============================
//interval tree stuff

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
void interval_tree::construct(vector<merged_interval *> D){
	if (D.size()>0){
		int N 		= D.size();
		int center 	= N / 2;
		current 	= D[center];
		ID 			= D[center]->id;
		//printf("%d,%d, %d,%d\n", ID, current->start, N, center);
		if ((center - 1) > 0){
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
	//printf("%d, %d\n", ID , N);
	ARRAY[ID] 	= *current;
	if (left!=NULL){
		left->insert_into_array(ARRAY, N);
	}
	if (right != NULL){
		right->insert_into_array(ARRAY, N);
	}
}




map<string, vector<merged_interval*> >  load_intervals(string FILE){

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
				G[chrom].push_back(interval(chrom, stoi(start), stoi(stop) ));
			
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
				intervals[chrom]->insert(double(i), coverage, strand);
			}		
		}

	}else{
		cout<<"Coudn't open: "<<FILE<<endl;

	}
}


void write_out(string FILE, map<string, interval_tree *> A){
	typedef map<string, interval_tree *>::iterator it_type;
	int N;
	ofstream FHW;
	
	FHW.open(FILE);
	
	for(it_type c = A.begin(); c != A.end(); c++) {
		N 	= A[c->first]->getDepth();
		merged_interval array[N];
		A[c->first]->insert_into_array(array,N);
		for (int i = 0; i < N;i++){
			for (int j = 0; j < array[i].intervals.size(); j++ ){
				if (array[i].intervals[j].forward_x.size()>0 and array[i].intervals[j].reverse_x.size()>0){
					FHW<<"#"<<c->first<<","<<to_string(array[i].intervals[j].start)<<","<<to_string(array[i].intervals[j].stop)<<endl;
					FHW<<"~forward"<<endl;
					for (int k = 0; k < array[i].intervals[j].forward_x.size(); k++ ){
						FHW<<to_string(int(array[i].intervals[j].forward_x[k]))<<","<<to_string(array[i].intervals[j].forward_y[k])<<endl;
					}
					FHW<<"~reverse"<<endl;
					for (int k = 0; k < array[i].intervals[j].reverse_x.size(); k++ ){
						FHW<<to_string(int(array[i].intervals[j].reverse_x[k]))<<","<<to_string(array[i].intervals[j].reverse_y[k])<<endl;
					}
				}
			}
		}
	}
}




















