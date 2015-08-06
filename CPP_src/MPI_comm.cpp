#include <mpi.h>
#include "MPI_comm.h"
#include "load.h"
#include <vector>
#include "across_segments.h"
using namespace std;

vector<segment *> slice_segments(vector<segment *> segments, int rank, int nprocs){

	int count 	= segments.size()/ nprocs;
	int start 	= rank * count;
    int stop 	= min(start + count, int(segments.size()));
    if (rank==(nprocs-1)){
    	stop 	= int(segments.size());
    }    
    vector<segment *> new_segs(segments.begin() + start, segments.begin()+stop);
    return new_segs;
}
int get_all_segs_size(vector<segment *> segments, int rank, int nprocs ){
	int SUM 	= segments.size();
	MPI_Status status;
	int sr;
	for (int job=1; job<nprocs; job++){
		MPI_Recv(&sr, 1, MPI_INT, job, 1, MPI_COMM_WORLD,&status);
		SUM+=sr;
	}
	return SUM;

}
void send_seg_size(vector<segment *> segments){
	int SUM 	= segments.size();
	MPI_Send(&SUM, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);	
}

map<int,int> get_all_bidir_sizes(vector<simple_c> fits, int nprocs ){
	MPI_Status status;
	int sr;
	map<int, int> table;
	for (int job=1; job<nprocs; job++){
		MPI_Recv(&sr, 1, MPI_INT, job, 1, MPI_COMM_WORLD,&status);
		table[job] 	= sr;
	}
	return table;
}

void send_bidir_size(vector<simple_c> fits){
	int SUM 	= fits.size();
	MPI_Send(&SUM, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);	
}



// void commit_MPI_simple_c_struct(MPI_Datatype mystruct){
// 	int blocklens[8]={1,1,1, 1 ,3,3,3,2};
// 	MPI_Datatype old_types[8]={MPI_INT, MPI_INT, MPI_DOUBLE, MPI_DOUBLE,
// 								MPI_DOUBLE, MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE};
// 	MPI_Aint indices[8] = {0, sizeof(int), sizeof(double), sizeof(double),
// 							sizeof(double), sizeof(double), sizeof(double), sizeof(double)  };
// 	MPI_Type_struct( 2, blocklens, indices, old_types, &mystruct );
// 	MPI_Type_commit( &mystruct );
	
// }
vector<simple_c> gather_all_simple_c_fits(vector<simple_c> root_simple_fits, 
	map<int,int> bidir_table, int N, int nproc ){ //for root, N expected number
	int count 	= N/nproc;
	vector<simple_c> all_fits ;
	for (int j =0; j < root_simple_fits.size(); j++){
		all_fits.push_back(root_simple_fits[j]);
	}
	MPI_Datatype mystruct;
	//int blocklens[8]={1,1,1, 1 ,3,3,3,2};
	int blocklens[4]={1, 1,4,9};
	
	// MPI_Datatype old_types[8]={MPI_INT, MPI_INT, MPI_DOUBLE, MPI_DOUBLE,
	// 							MPI_DOUBLE, MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE};
	MPI_Datatype old_types[4] 	={MPI_DOUBLE, MPI_DOUBLE,MPI_INT, MPI_DOUBLE}; 
	MPI_Aint displacements[4];
	MPI_Aint intex, doublex;

	MPI_Type_extent(MPI_INT, &intex);

	MPI_Type_extent(MPI_DOUBLE, &doublex);



	displacements[0] = (MPI_Aint) 0;
	displacements[1] = doublex;
	displacements[2] = doublex+doublex;
	displacements[3] = doublex+doublex+intex*4 ;
	
	MPI_Type_struct( 4, blocklens, displacements, old_types, &mystruct );
	MPI_Type_commit( &mystruct );
	
	MPI_Status status;
	simple_c sr;
	typedef map<int,int>::iterator it_type;
	map<int, map<int, bidir_preds> > G; //map refering to segment ID and specific bidir ID
	map<int, bidir_preds> A;

	typedef map<int, map<int, bidir_preds> >::iterator it_type_2;
	typedef map<int, bidir_preds>::iterator it_type_3;


	int segment_id, bidir_id, Complexity, pred_bidirs_in_merged;
	for (it_type cc=bidir_table.begin(); cc!=bidir_table.end(); cc++){
		//cc->first is the job that we are expect it from
		for (int j = 0; j < cc->second ; j++){
			MPI_Recv(&sr, 1, mystruct, cc->first, j, MPI_COMM_WORLD,&status);
			
			segment_id=sr.IDS[0]+cc->first*count , bidir_id=sr.IDS[3], Complexity=sr.IDS[2], pred_bidirs_in_merged=sr.IDS[1];
			if (G.find(segment_id) == G.end()  ){
				G[segment_id] 	= A;
			}
			if (G[segment_id].find(bidir_id) == G[segment_id].end() ){
				G[segment_id][bidir_id] 	= bidir_preds(sr.noise_ll);
			}
			G[segment_id][bidir_id].insert_component(Complexity, sr, sr.ll);
		}
	}
	for (it_type_2 s = G.begin(); s!=G.end(); s++){
		printf("Segment: %d\n", s->first );
		for (it_type_3 b = s->second.begin(); b!=s->second.end(); b++ ){
			printf("Bidir: %d\n", b->first);
		}
	}

	

	return all_fits;
}

void send_all_simple_c_fits(vector<simple_c> root_simple_fits ){ //slaves
	//need to setup MPI derived datatype for simple_c struct
	MPI_Datatype mystruct;
	int blocklens[4]={1, 1,4,9};
	MPI_Datatype old_types[4] = {MPI_DOUBLE, MPI_DOUBLE,MPI_INT, MPI_DOUBLE}; 
	MPI_Aint displacements[4];
	MPI_Aint intex, doublex;

	MPI_Type_extent(MPI_INT, &intex);

	MPI_Type_extent(MPI_DOUBLE, &doublex);



	displacements[0] = (MPI_Aint) 0;
	displacements[1] = doublex;
	displacements[2] = doublex+doublex;
	displacements[3] = doublex+doublex+intex*4 ;
	
	

	MPI_Type_struct( 4, blocklens, displacements, old_types, &mystruct );
	MPI_Type_commit( &mystruct );
	
	for (int j = 0; j < root_simple_fits.size(); j++){
		MPI_Send(&root_simple_fits[j], 1, mystruct, 0, j, MPI_COMM_WORLD);
	}




}

void gather_all_segments(){ //for root

}

void send_all_segments(){ //slaves

}









