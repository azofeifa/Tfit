#ifndef load_H
#define load_H
#include <string>
#include <map>
#include "pdfs.h"

using namespace std;

class annotation{
public:
	int start, stop,N, step;
	double total;
	string name, strand, chrom;
	double ** D;
	annotation * next; 
	double density;
	DGU * rv_dgu_best;
	void add(double, double, double, int);
	annotation(string, string, string, int, int, int);
	annotation();
	~annotation();
};
class annotation_cluster{
public:
	int start, stop;
	annotation * cluster;
	annotation * root;
	
	annotation_cluster* next;
	annotation_cluster* current;
	annotation_cluster(annotation *);
	~annotation_cluster();
	void add_interval(annotation *);
	void find(double, double , double , int);

};

class annotations{
public:
	map<string, annotation_cluster *> collections;
	double N;
	annotations();
	~annotations(void);
	void add_annotation(string, string, string, int, int, int);
	void add_coverage(string, double, double, double, int);
};


int read_bedgraph(string, annotations *, int);
int read_annotation(string, int, annotations *,int);
void standardize(annotations *);
#endif