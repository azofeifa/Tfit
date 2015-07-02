#include "load.h"
#include <string>
#include <iostream>
#include <fstream>
#include "split.h"
#include "pdfs.h"
using namespace std;
annotation::annotation(string n, string s, string c, int st, int sp, int r){
	name 		= n;
	strand 		= s;
	chrom  		= c;
	start 		= st;
	stop 		= sp;
	step 		= r;
	N 			= (sp -st) / r;
	next 		= NULL;
	total 		= 0;
	density		= 0;
	D 			= NULL;
	rv_dgu_best = NULL;
}
void annotation::add(double start, double stop, double cov, int s){
	if (D==NULL){
		D  		= new double*[N];
		int st 	= start + step;
		for (int i = 0; i < N; i++){
			D[i] 		= new double[3];
			D[i][0] 	= (double) st;
			D[i][1] 	= 0;
			D[i][2] 	= 0;
			st 			= st + step;
		}
	}
	int i = 0;
	for (int j = start; j < stop+1; j++){
		while (i < N and D[i][0] < j){
			i+=1;
		}
		if (i<N){
			D[i][s]+=cov;
			total+=cov;	
		}
	}
}
annotation_cluster::annotation_cluster(annotation * C){
	start 	= C->start;
	stop 	= C->stop;
	cluster = C;
	next 	= NULL;
	current = this;
	root 	= C;
}
annotation_cluster::~annotation_cluster(){
	annotation * current 	= root;
	annotation * next;
	while(current!=NULL){
		next 	= current->next;
		delete(current);
		current = next;
	}

}

void annotation_cluster::add_interval(annotation * C){
	if (C->start < start)	{
		start 	= C->start;
	}
	if (C->stop > stop){
		stop 	= C->stop;
	}
	cluster->next 	= C;
	cluster 		= cluster->next;

}

void annotation_cluster::find(double st, double sp, double cov,int  s){
	annotation * sec_current 	= NULL;
	while(current!=NULL and current->stop  < st){
		current 	= current->next;
	}
	if (current and st<=current->stop and sp >= current->start ){
		sec_current 			= current->root;
	}
	if (sec_current!= NULL){
		while (sec_current!=NULL) {
			if (st<=sec_current->stop and sp >= sec_current->start){
				sec_current->add(st, sp, cov, s);
			}
			sec_current 	= sec_current->next; 
		}
	}
}

annotation::~annotation(){
	if (total>0){
		for (int i =0; i < N; i++){
			delete D[i];
		}
		delete D;
	}

}
annotation::annotation(){
	total 	= 0;
	N 		= 0;

}

annotations::annotations(){
	N=0;
};

annotations::~annotations(){
	cout<<"cleaning up/exiting..."<<endl;
	typedef map<string, annotation_cluster *>::iterator it_type;
	annotation_cluster * cluster;
	annotation_cluster * c_next; 
		
	for(it_type iterator 	= collections.begin(); iterator != collections.end(); iterator++) {
		cluster 			= iterator->second;
		c_next 				= cluster;
		while (cluster!=NULL){
			c_next 					= cluster->next;
			delete(cluster);
			cluster 				= c_next;
		}
	}
	
}

int max(int x, int y){
	if (x>y){
		return x;
	}
	return y;
}
int min(int x, int y){
	if (x<y){
		return x;
	}
	return y;
}


void annotations::add_annotation(string name, string strand, string chrom, int start, int stop, int r){
	annotation * insertee 			= new annotation(name, strand, chrom, start, stop, r);
	N+=1;
	if (collections.find(chrom) == collections.end()){
		collections[chrom] 	 = new annotation_cluster( insertee );
	}else{
		annotation_cluster * current 	= collections[chrom];
		annotation_cluster * prev 		= NULL;
		while(current!=NULL and current->stop < start){
			prev 		= current;
			current 	= current->next;
		}
		annotation_cluster * AC 	= new annotation_cluster(insertee);
		if (current==NULL){
			prev->next 	= AC;
		}else if(stop < current->start){
			AC->next 	= current;
			if (prev==NULL){
				collections[chrom] 	= AC;	
			}else{
				prev->next 	= AC;
			}
		}
		else{
			annotation * c_cluster;
			int u_st 	= min(start, current->start);
			int u_sp 	= max(stop, current->stop);
			while(current!=NULL and current->start < stop){
				u_st 	= min(u_st, current->start);
				u_sp 	= max(u_sp, current->stop);
				c_cluster 	= current->root;
				while (c_cluster!=NULL){
					AC->cluster->next 	= c_cluster;
					AC->cluster 		= AC->cluster->next;
					c_cluster 			= c_cluster->next;
				}
				current 	= current->next;
			}
			AC->start 	= u_st;
			AC->stop 	= u_sp;

			AC->next 	= current;
			if (prev==NULL){
				collections[chrom] 	= AC;	
			}else{
				prev->next 	= AC;
			}

		}
	}
}

void annotations::add_coverage(string chrom, double st, double sp, double cov, int s){
	if (collections.find(chrom) != collections.end()){
		collections[chrom]->find(st, sp, cov, s);
	}
}

int read_annotation(string FILE, int r, annotations * A, int pad){
	ifstream FH(FILE);
	int N 	= 0;
	if (FH){
		string line, chrom, start, stop, strand, name; 
		vector<string> lineArray;
		while (getline(FH, line)){
			if ("#"!=line.substr(0,1)){
				lineArray 	= splitter(line, "\t");
				name 		= lineArray[1];
				chrom 		= lineArray[2];
				strand 		= lineArray[3];
				start 		= lineArray[4];
				stop 		= lineArray[5];
				A->add_annotation(name, strand, chrom, stoi(start)-pad, stoi(stop)+pad, r);
				N++;
			}
		}
	}else{
		cout<<"could not open: "<<FILE<<" "<<endl;
		return 0;
	}
	return 1;
}
void set_roots(annotations * A){
	map<string, annotation_cluster *> C 	= A->collections;
	typedef map<string, annotation_cluster *>::iterator it_type;
	annotation_cluster * cluster;
	for(it_type iterator 	= C.begin(); iterator != C.end(); iterator++) {
		cluster 			= iterator->second;
		cluster->current 	= cluster;
	}
}

int read_bedgraph(string FILE, annotations * A, int s){
	set_roots(A);
	ifstream FH(FILE);
	if (FH){
		string line, chrom, start, stop, cov;
		vector<string> lineArray;
		while (getline(FH, line)){
			lineArray 	= splitter(line, "\t");
			chrom 		= lineArray[0];
			start 		= lineArray[1];
			stop 		= lineArray[2];
			cov 		= lineArray[3];
			A->add_coverage(chrom, stod(start), stod(stop), stod(cov), s );
		}
	}else{
		cout<<"could not open: " << FILE<<endl;
		return 0;
	}	
	map<string, annotation_cluster *> C 	= A->collections;
	typedef map<string, annotation_cluster *>::iterator it_type;
	double KB  = 1000;
	int start;		
	annotation_cluster * cluster;
	annotation * current;
	return 1;
}

void standardize(annotations * A){
	map<string, annotation_cluster *> C 	= A->collections;
	typedef map<string, annotation_cluster *>::iterator it_type;
	double KB  = 1000;
	int start;		
	annotation_cluster * cluster;
	annotation * current;
	for(it_type iterator 	= C.begin(); iterator != C.end(); iterator++) {
		cluster 			= iterator->second;
		while (cluster!=NULL){
			annotation * current 	= cluster->root;
			while (current!=NULL){
				if (current->total > 0){
					current->density 	= current->total / (current->stop - current->start);
					//get rid of zero entries
					int NN 	= 0;
					for (int i =0; i < current->N;i++){
						if (current->D[i][1] or current->D[i][2] ){
							NN+=1;
						}
					}
					double ** DD;
					DD 	= new double*[NN];
					int j 	= 0;
					for (int i = 0; i < current->N; i++ ){
						if (current->D[i][1] or current->D[i][2] ){
							DD[j] 		= new double[3];
							DD[j][0] 	= current->D[i][0];
							DD[j][1] 	= 0;
							DD[j][2] 	= 0;
							if (current->D[i][1]){
								DD[j][1]= current->D[i][1];
							}
							if (current->D[i][2]){
								DD[j][2]= current->D[i][2];
							}
							j++;
						} 
					}
					current->D 	= DD;
					current->N 	= NN;
					//center to minimal value and standardize
					int MIN 				= current->D[0][0];
					for (int i = 0; i < current->N; i++){
						current->D[i][0]	= (current->D[i][0]);
					}
				}
				current=current->next;
			}
			cluster 	= cluster->next;
		}
	}
}

