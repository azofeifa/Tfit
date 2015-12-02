#ifndef MPI_comm_H
#define MPI_comm_H
#include "load.h"
#include "across_segments.h"
#include <vector>
#include <fstream>
#include <iostream>

#include <map>
#include "read_in_parameters.h"
namespace MPI_comm {

 
vector<segment *> slice_segments(vector<segment *>, int , int );

int gather_all_bidir_predicitions(vector<segment *> ,
vector<segment *>, int, int,string, string, int, params *, int );


map<string, vector<segment *> > send_out_single_fit_assignments(vector<segment *> , int, int);


int get_job_ID(string,string,int, int);


map<int, map<int, vector<simple_c_free_mode>  > >  gather_all_simple_c_free_mode(vector<map<int, vector<simple_c_free_mode> >>  , 
	int  , int );

void wait_on_root(int, int);

}
#endif
