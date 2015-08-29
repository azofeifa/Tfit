#ifndef MPI_comm_H
#define MPI_comm_H
#include "load.h"
#include "across_segments.h"
#include <vector>
#include <map>
struct rsimple_c{
public:
	int st_sp[5]; //start and stop of the bidir segment
	char chrom[5];
	//need some important IDS
	//1: the segment from all_segments that this bidir segment came from, mainly for the chromosome ID
	//2: we need to know what bidir this component refers
	//3: need to konw what the complexity is so that we can match it up with other components that have the same stats
	//####
	//below is just all the parameter estimates for the component, EMG and forward/reverse uniforms
	double ps[14]; //0->mu, 1->si, 2->l, 3->w, 
				   //4->pi, 5->fw, 6->rw, 7->fb, 8->ra, 9->fpi, 10->rpi, 11->NN
	rsimple_c();
};
vector<segment *> slice_segments(vector<segment *>, int , int );
int get_all_segs_size(vector<segment *>, int , int);
void send_seg_size(vector<segment *>);
map<int,int> get_all_bidir_sizes(vector<simple_c>, int);
void send_bidir_size(vector<simple_c>);
map<int, map<int, bidir_preds> > gather_all_simple_c_fits(vector<simple_c>, map<int,int>, int , int );
void send_all_simple_c_fits(vector<simple_c>);
map<string , vector<vector<double> > > gather_all_bidir_predicitions(vector<segment *> ,
vector<segment *>, int, int,string );
map<string, map<int, vector<rsimple_c> > > gather_all_simple_c_fits(vector<segment *>, vector<simple_c>, int , int);

#endif
