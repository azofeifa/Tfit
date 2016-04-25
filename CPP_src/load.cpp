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
#ifdef USING_ICC
#include <mathimf.h>
#else
#include <math.h>   
#endif
#include <limits>
#include <iostream>
#include <algorithm>
#include <unistd.h>
#include <random>
#include <exception>

#include <stdio.h>
#include <time.h>
#ifdef USING_ICC
#include <mathimf.h>
#include <aligned_new>
#else
#include <math.h>   
#endif

using namespace std;
//========================================================================
//The very very important segment class

segment::segment(string chr, int st, int sp){
	chrom	= chr;
	start	= st;
	stop	= sp;
	N 		= 0;
	fN 		= 0;
	rN 		= 0;
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
	fN 		= 0;
	rN 		= 0;
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
	fN 		= 0;
	rN 		= 0;
	minX=st, maxX=sp;
	counts 	= 0;
	XN 		= 0;
	ID 		= i;
	strand 	= STR;
	chrom_ID= 0;

}

segment::segment(){
	N 		= 0;
	fN 		= 0;
	rN 		= 0;
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
	fN = 0, rN = 0;
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
			fN+=forward[i][1];
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
			rN+=reverse[i][1];
		}
	}
	//===================
	//scale data down for numerical stability
	if (scale){
		for (int i = 0; i < BINS; i ++ ){

			X[0][i] 	= (X[0][i]-minX)/scale;
			// X[1][i]/=delta;
			// X[2][i]/=delta;
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

//================================================================================================
//interval tree code

node::node(){};

node::node(vector<segment * > segments ){
	center 	= (double(segments[0]->start)  + double(segments[segments.size()-1]->stop)) / 2.;
	vector<segment * > Left;
	vector<segment * > Right;
	left=NULL, right=NULL;
	for (int i = 0 ; i < segments.size(); i++){
		if (segments[i]->stop < center){
			Left.push_back(segments[i]);
		}
		else if (segments[i]->start > center){
			Right.push_back(segments[i]);
		}
		else{
			current.push_back(segments[i]);
		}
	}
	if (Left.size() > 0){
		left 	= new node(Left);
	}
	if (Right.size() > 0){
		right 	= new node(Right);
	}
}
void node::insert_coverage(vector<double> x, int s){
	for (int i = 0 ; i < current.size(); i++){
		if (x[0] > current[i]->start and  x[0] < current[i]->stop  ){
			if (s==1){
				current[i]->forward.push_back(x);
			}else{
				current[i]->reverse.push_back(x);	
			}
		}
	}	

	if (x[0] >= center and right != NULL ){
		right->insert_coverage(x, s);
	}
	if (x[0] <= center and left !=NULL){
		left->insert_coverage(x,  s);
	}
}
void node::searchInterval(int start, int stop, vector<int>& finds ){
	for (int i = 0 ; i < current.size(); i++){
		if (stop > current[i]->start and  start < current[i]->stop  ){
			finds.push_back(1);
		}
	}	
	if (start >= center and right != NULL ){
		right->searchInterval(start, stop, finds);
	}
	if (stop <= center and left !=NULL){
		left->searchInterval(start, stop, finds);
	}	
}

void node::retrieve_nodes(vector<segment*> & saves){
	for (int i = 0; i < current.size(); i++){
		saves.push_back(current[i]);
	}
	if (right!= NULL){
		right->retrieve_nodes(saves);
	}
	if (left != NULL){
		left->retrieve_nodes(saves);		
	}
}

//================================================================================================
//a class from loading the K_models formmated file 
segment_fits::segment_fits(){} //empty con
segment_fits::segment_fits(string c, int st, int sp,
 	double n_pos, double n_neg, string id){

	chrom=c, start=st, stop=sp, TSS=0;
	N 	= n_pos + n_neg, N_pos = n_pos, N_neg = n_neg;
	ID 	= id;
}
void segment_fits::get_model(double ms_pen){
	typedef map<int, double>::iterator it_type;
	int arg;
	double MIN=INF;
	double score;
	
	for (it_type m 	= M.begin(); m != M.end(); m++){
		if (m->first > 0){
			score 	= -2*m->second + log(N)*(ms_pen*m->first);
		}else{
			score 	= -2*m->second + log(N) ;		
		}
		if (score < MIN){
			MIN 	= score;
			arg 	= m->first;
		}
	}
	model 	= arg;
}
string segment_fits::write (){
	string line 				= "";
	if (model > 0){
		string forward_N=to_string(N_pos), reverse_N = to_string(N_neg);
		vector<string> params 		= split_by_tab(parameters[model], "");
		for (int i = 0 ; i < model; i++){
			vector<string> S 	=  split_by_comma(params[0], "");
			double mu 	= stod(split_by_comma(params[0], "")[i]);
			double std 	= stod(split_by_comma(params[1], "")[i]);
			double lam 	= stod(split_by_comma(params[2], "")[i]);
			double w 	= stod(split_by_comma(split_by_bar(params[5], "")[i] , "" )[0] );

			int start 	= mu - (std + lam), stop = mu + (std+lam);
			if (std  < 50000 and lam < 20000 and w > 0.05 ){
				line+=chrom+"\t" + to_string(start) + "\t" + to_string(stop)+"\t";
				line+=ID+"|";
				string ps 	= "";
				for (int m = 0; m < params.size(); m++){
					vector<string> ks;
					if (m!=5){
						ks 	= split_by_comma(params[m], "");
					}else{
						ks 	= split_by_bar(params[m], "");	
					}
					if (m + 1 < params.size()){
						ps+=ks[i]+ ",";
					}else{
						ps+=ks[i];
					}
				}
				line+=ps+"\n";
			}
		}
	}
	return line;
}

//================================================================================================
//merge segments from loading_intervals
vector<segment *> merge_segments(vector<segment *> segments, map<int, string>  IDS_first, map<int, string> & IDS, int & T){
	vector<segment *> new_segments;
	//bubble sort
	bool changed 	= true;
	while (changed){
		changed = false;
		for (int i = 1 ; i < segments.size(); i++){
			if (segments[i-1]->start > segments[i]->start){
				changed 		= true;
				segment * copy 	= segments[i-1];
				segments[i-1] 	= segments[i];
				segments[i] 	= copy;
			}
		}
	}

	int j = 0, N = segments.size(), i =0;
	while (j < N){
		string ID 		= "";
		segment * S 	= new segment(segments[i]->chrom, 
			segments[j]->start,segments[j]->stop, T, segments[j]->strand );
		while (j < N and segments[j]->start < S->stop and segments[j]->stop > S->start ){
			S->start 	= min(S->start,segments[j]->start);
			S->stop 	= max(S->stop, segments[j]->stop);
			ID+=IDS_first[segments[j]->ID] + ",";
			S->counts+=1;
			j++;
		}
		new_segments.push_back(S);
		IDS[T] 			= ID.substr(0, ID.size()-1);
		T+=1;
	}

	return new_segments;
}
bool check_ID_name(string & INFO){
	bool PASSED 	= true;
	string change 	= "::";
	for (int i = 0; i < INFO.size(); i++){
		if (INFO.substr(i,1)=="|"){
			PASSED 		= false;
			INFO.replace(i, 1, change);
		}
	}
	return PASSED;
}

//================================================================================================
//LOADING from file functions...need to clean this up...





vector<segment*> load::load_bedgraphs_total(string forward_strand, 
	string reverse_strand, string joint_bedgraph, int BINS, double scale, string spec_chrom, map<string, int>& chromosomes
	, map<int, string>& ID_to_chrom){
	bool FOUND 	= false;
	if (spec_chrom=="all"){
		FOUND 	= true;
	}
	map<string, segment*> 	G;
	vector<segment*> segments;
	vector<string> FILES;
	if (forward_strand.empty() and reverse_strand.empty()){
		FILES 	= {joint_bedgraph};
	}else if (not forward_strand.empty() and not reverse_strand.empty()){
		FILES 	= {forward_strand, reverse_strand};
	}
	
	string line, chrom;
	int start, stop;
	double coverage;
	vector<string> lineArray;
	string prevChrom="";
	segment * S =NULL;
	bool INSERT 	= false;
	bool EXIT 		= false;
	int line_number = 0;
	for (int u = 0 ; u < FILES.size(); u++){
		ifstream FH(FILES[u]) ;
		if (not FH ){
			printf("couln't open FILE %s\n", FILES[u].c_str());
		}
		if (EXIT){
			break;
		}
		while (getline(FH, line)){
			lineArray=splitter(line, "\t");
			if (lineArray.size()!=4){
				EXIT 	= true;
				printf("\nLine number %d  in file %s was not formatted properly\nPlease see manual\n",line_number, FILES[u].c_str() );
			}
			line_number++;
			chrom=lineArray[0], start=stoi(lineArray[1]), stop=stoi(lineArray[2]), coverage=(stof(lineArray[3]));
			if (chrom != prevChrom and (chrom==spec_chrom or spec_chrom=="all")  )  {
				FOUND 		= true;
				if (chrom.size() < 6 and u==0){
					G[chrom] 	= new segment(chrom, start, stop );
					INSERT 		= true;
					FOUND 		= true;
				}else if(chrom.size() > 6){
					INSERT 		= false;
				}
			}
			if (FOUND and chrom!= spec_chrom and spec_chrom!= "all"){
				break;
			}
			if (INSERT){
				if (u==0){
					if (coverage > 0){
						G[chrom]->add2(1, double(stop+start)/2., abs(coverage));
					}else{
						G[chrom]->add2(-1, double(stop+start)/2., abs(coverage));						
					}
				}else{
					G[chrom]->add2(-1, double(stop+start)/2., abs(coverage));	
				}
			}
			prevChrom=chrom;

		}
	}
	if (not EXIT){
		int c =1;
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
	}
	if (not FOUND){
		segments.clear();
		printf("couldn't find chromosome %s in bedgraph files\n", spec_chrom.c_str());
	}
	return segments;
}
vector<segment*> load::load_intervals_of_interest(string FILE, map<int, string>&  IDS, 
	params * P, bool center){
	ifstream FH(FILE);

	string spec_chrom 	= P->p["-chr"];
	int pad 			= stoi(P->p["-pad"])+1;
	int merge 			= stof(P->p["-merge"]);

	vector<segment *> G;
	int ct 	= 1;
	map<string, vector<segment * > > GS;
	map<int, string> IDS_first;
	int T 	= 0;
	bool EXIT 		= false;
	if (FH){
		string line, chrom;
		int start, stop;
		int 	i = 0;
		vector<string> lineArray;
		string strand; 
		bool PASSED 	= true;

		while(getline(FH, line)){
			lineArray=splitter(line, "\t");
			if (lineArray[0].substr(0,1)!="#" and lineArray.size()>2){
				if (lineArray.size() > 3){
					if (not check_ID_name(lineArray[3]) and PASSED ){
						PASSED 			= false;
						printf("\ninterval id in line: %s, contains a | symbol changing to :: -> %s\n",line.c_str(), lineArray[3].c_str() );
						printf("Will continue to change other occurrences....\n");

					}
					IDS_first[i] 		= lineArray[3];
				}else{
					IDS_first[i] 		= "Entry_" + to_string(i+1);	
				}
				if (lineArray.size() > 4){
					strand 		= lineArray[4];
				}else{
					strand 		= ".";
				}
				try{
					if (not center){
						chrom=lineArray[0], start=max(stoi(lineArray[1])-pad, 0), stop=stoi(lineArray[2]) + pad;
					}else{
						int x 	= 	((stoi(lineArray[1]) + stoi(lineArray[2])))/2.;
						start 		= max(x - pad, 0) , stop 	= x + pad;
						chrom=lineArray[0];
					}
				}
				catch(exception& e){
					printf("\n\nIssue with file %s at line %d\nPlease consult manual on file format\n\n",FILE.c_str(), i );
					EXIT=true;
					GS.clear();
					break;
				}
				if (start < stop){
					if (spec_chrom=="all" or spec_chrom==chrom){

						segment * S 	= new segment(chrom, start, stop,i,strand);
						GS[S->chrom].push_back(S);
					}
					i++;
				}
			}
		}
	}else{
		printf("couldn't open %s for reading\n", FILE.c_str() );
		EXIT 	= true;
	}
	if (not EXIT){
		typedef map<string, vector<segment * > >::iterator it_type;
		if (not merge){
			IDS 	= IDS_first;
		}
		for (it_type c 	= GS.begin(); c!=GS.end(); c++){
			vector<segment *> m_segs;
			if (merge){
				m_segs 	= merge_segments(c->second, IDS_first, IDS, T);
			}else{
				m_segs 	= c->second;
			}
			for (int i = 0 ; i < m_segs.size(); i++){
				G.push_back(m_segs[i]);
			}
		}
	}else{
		G.clear();
	}
	return G;
}

vector<segment* > load::insert_bedgraph_to_segment_joint(map<string, vector<segment *> > A , 
	string forward, string reverse, string joint, int rank ){
	
	
	
	map<string, node> NT;
	typedef map<string, vector<segment *> >::iterator it_type_5;
	for(it_type_5 c = A.begin(); c != A.end(); c++) {
		NT[c->first] 	= node(c->second);
	}
	int start, stop, N, j;
	double coverage;
	N 	= 0,j 	= 0;
	int strand;
	int o_st, o_sp;
	vector<string> lineArray;
	string chrom, prevchrom, line;
	vector<segment *> segments;
	double center;
	vector<string> FILES;
	if (forward.empty() and reverse.empty()){
		FILES 	= {joint};
	}else if (not forward.empty() and not reverse.empty()) {
		FILES 	= {forward, reverse};
	}
	string FILE;
	for (int i =0; i < FILES.size(); i++){
		FILE=FILES[i];
		ifstream FH(FILE);
		if (FH){
			prevchrom="";
			while (getline(FH, line)){
				lineArray 	= splitter2(line, "\t");
				if (lineArray.size()==4){
					chrom 		= lineArray[0];
					start=stoi(lineArray[1]),stop=stoi(lineArray[2]), coverage = stod(lineArray[3]);
					if (coverage > 0 and i == 0){
						strand 	= 1;
					}else if (coverage < 0 or i==1){
						strand 	= -1;
					}
					center 	= (stop + start) /2.;
					if (NT.find(chrom)!=NT.end()){
						vector<double> x(2);
						x[0]=center, x[1] = abs(coverage);
						NT[chrom].insert_coverage(x, strand);
						
					}
				}
				else{
					printf("\n***error in line: %s, not bedgraph formatted\n", line.c_str() );
					segments.clear();
					return segments;
				}
			}
			FH.close();
			
		}else{
			cout<<"could not open forward bedgraph file: "<<FILE<<endl;
			segments.clear();
			return segments;
		}
	}
	//now we want to get all the intervals and make a vector<segment *> again...
	vector<segment *>NS;
	typedef map<string, node>::iterator it_type_6;
	for (it_type_6 c = NT.begin(); c!=NT.end(); c++){
		c->second.retrieve_nodes(NS);
	}

	return NS;
}
vector<segment_fits *> load::load_K_models_out(string FILE){
	ifstream FH(FILE);
	string line;
	segment_fits * S = NULL;
	int complexity;
	double ll;
	string chrom;
	int start,stop;
	vector<segment_fits *> segment_fits_all;	
	if (FH){
		while (getline(FH, line)){
			if (line.substr(0,1)!= "#" and line.substr(0,1)==">"){
				if (S!=NULL){
					segment_fits_all.push_back(S);
				}
				line 							= line.substr(1,line.size()-1);

				vector<string> bar_split 		= split_by_bar(line, "");
				vector<string> comma_split 		= split_by_comma(bar_split[2], "");
				vector<string> colon_split 		= split_by_colon(bar_split[1], "");
				vector<string> dash_split 		= split_by_dash(colon_split[1], "");
				chrom = colon_split[0], start 	= stoi(dash_split[0]) ,stop = stoi(dash_split[1]) ;
				S 		= new segment_fits(chrom,start,
				 stop, stod(comma_split[0]), stod(comma_split[1]), bar_split[0] );

			}else if (line.substr(0,1)=="~" and S!=NULL){
				line 	= line.substr(1,line.size()-1);
				vector<string> tab_split 		= split_by_tab(line, "");
				vector<string> comma_split 		= split_by_comma(tab_split[0], "");
				complexity 	= stoi(comma_split[0]), ll 	= stod(comma_split[1]);
				S->M[complexity]=ll;
				string parameters 	= "";
				for (int i = 1; i < tab_split.size(); i++){
					if (i+1 < tab_split.size()){
						parameters+=tab_split[i]+"\t";
					}else{
						parameters+=tab_split[i];
					}
				}
				S->parameters[complexity] 	= parameters;
				
			}
		}
		if (S!=NULL){
			segment_fits_all.push_back(S);		
		}
	}else{
		printf("couldn't open %s...strange error\n", FILE.c_str() );
	}
	return segment_fits_all;
}


//================================================================================================
//WRITE out to file functions


void load::write_out_bidirs(map<string , vector<vector<double> > > G, string out_dir, 
	string job_name,int job_ID, params * P, int noise){
	typedef map<string , vector<vector<double> > >::iterator it_type;
	ofstream FHW;
	if (not noise){
		FHW.open(out_dir+ job_name+ "-" + to_string(job_ID)+ "_prelim_bidir_hits.bed");
	}else{
		FHW.open(out_dir+ job_name+ "-" + to_string(job_ID)+ "_non_div_txn.bed");	
	}
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


void load::write_out_models_from_free_mode(map<int, map<int, vector<simple_c_free_mode>  > > G, 
	params * P, int job_ID,map<int, string> IDS, int noise, string & file_name){

	//========================================================================================
	//write out each model parameter estimates
	double scale 	= stof(P->p["-ns"]);
	double penality = stof(P->p["-ms_pen"]) ;
	string out_dir 	= P->p["-o"];
	ofstream FHW;
	if (noise==0){
		file_name 	= out_dir+  P->p["-N"] + "-" + to_string(job_ID)+  "_K_models_MLE.tsv";
		FHW.open(file_name);
	}else{
		file_name 	= out_dir+  P->p["-N"] + "-" + to_string(job_ID)+  "_non_div_txn_K_models_MLE.tsv";
		FHW.open(file_name);
	}
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
	FHW.flush();
}
void load::write_out_bidirectionals_ms_pen(vector<segment_fits*> fits, params * P, int job_ID, int noise ){
	ofstream FHW;
	if (noise){
		FHW.open(P->p["-o"]+  P->p["-N"] + "-" + to_string(job_ID)+  "_non_div_txn_divergent_classifications.bed");
	}else{
		FHW.open(P->p["-o"]+  P->p["-N"] + "-" + to_string(job_ID)+  "_divergent_classifications.bed");	
	}
	FHW<<P->get_header(0);
	double penality 	= stod(P->p["-ms_pen"]);
	for (int i = 0; i < fits.size(); i++){
		fits[i]->get_model(penality);
		FHW<<fits[i]->write();
	}
	FHW.flush();
}

//================================================================================================
//misc.
void load::BIN(vector<segment*> segments, int BINS, double scale, bool erase){
	for (int i = 0 ; i < segments.size() ; i ++){
		if (segments[i]->forward.size() > 0 or segments[i]->reverse.size() > 0 ){
			segments[i]->bin(BINS, scale, erase);
		}
	}
}

void load::collect_all_tmp_files(string dir, string job_name, int nprocs, int job_ID){
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
		}
	}
}
void load::clear_segments(vector<segment *> segments){
	for (int i = 0; i < segments.size(); i++){
		if (segments[i]!=NULL){
			delete (segments[i]);
		}
	}
}
vector<segment_fits *> load::label_tss(string tss_file, vector<segment_fits *> query_fits ){
	vector<segment_fits *> new_fits;
	ifstream FH(tss_file);
	string line;
	map<string, vector<segment *> >G;
	map<string, node> T;
	string chrom,start, stop;
	typedef map<string, vector<segment *>>::iterator it_type;
	if (FH){	
		while (getline(FH, line)){
			vector<string> lineArray 	= split_by_tab(line, "");
			chrom 	= lineArray[0], start = lineArray[1], stop = lineArray[2];
			segment * S 	= new segment(chrom, stoi(start), stoi(stop));
			G[chrom].push_back(S);
		}
		//make G a node interval tree
		for (it_type c = G.begin(); c!=G.end(); c++){
			T[c->first] 	= node(c->second);
		}
		//label
		for (int i = 0 ; i < query_fits.size(); i++){
			vector<int> FINDS;
			if (T.find(query_fits[i]->chrom) != T.end()){
				T[query_fits[i]->chrom].searchInterval(query_fits[i]->start,query_fits[i]->stop, FINDS);
				if (!FINDS.empty() and query_fits[i]->M[1]!=nINF and query_fits[i]->M[1] > query_fits[i]->M[0]   ){
					query_fits[i]->TSS=1;
					new_fits.push_back(query_fits[i]);
				}
			}
			FINDS.clear();

		}
		
	}else{
		printf("couldn't open tss file %s... weird error\n", tss_file.c_str() );
	}
	return new_fits;	
}
















