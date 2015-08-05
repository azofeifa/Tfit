#ifndef MPI_comm_H
#define MPI_comm_H
#include "load.h"
#include "across_segments.h"
#include <vector>
#include <map>
vector<segment *> slice_segments(vector<segment *>, int , int );
int get_all_segs_size(vector<segment *>, int , int);
void send_seg_size(vector<segment *>);
map<int,int> get_all_bidir_sizes(vector<simple_c>, int);
void send_bidir_size(vector<simple_c>);
vector<simple_c> gather_all_simple_c_fits(vector<simple_c>, map<int,int>, int , int );
void send_all_simple_c_fits(vector<simple_c>);
#endif
