#include <mpi.h>
#include "MPI_comm.h"
#include "load.h"
#include <vector>
#include "across_segments.h"
#include <sys/types.h>
#include <dirent.h>
#include "split.h"
#include <stddef.h>
#include "read_in_parameters.h"
using namespace std;

vector<segment *> slice_segments(vector<segment *> segments, int rank, int nprocs){

	int count 	= segments.size()/ nprocs;
	if (count==0){
		count 	= 1;
	}
	int start 	= rank * count;
    int stop 	= min(start + count, int(segments.size()));

    if (rank==(nprocs-1)){
    	stop 	= int(segments.size());
    }
    if (start >= stop){
    	start 	= stop;
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

double send_density_val(double density, int rank, int nprocs){
	double D;
	MPI_Status status;
	if (rank==0){
		D 	= density;
		for (int j = 1; j < nprocs; j++){
			MPI_Send(&D, 1, MPI_DOUBLE,j,0, MPI_COMM_WORLD);
		}
	}else{
		MPI_Recv(&D, 1, MPI_DOUBLE,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE )	;
	}
	return D;


}

void send_bidir_size(vector<simple_c> fits){
	int SUM 	= fits.size();
	MPI_Send(&SUM, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);	
}



void send_all_simple_c_fits(vector<simple_c> root_simple_fits ){ //slaves
	//need to setup MPI derived datatype for simple_c struct
	MPI_Datatype mystruct;
	int blocklens[4]={1, 1,4,12};
	MPI_Datatype old_types[4] = {MPI_DOUBLE, MPI_DOUBLE,MPI_INT, MPI_DOUBLE}; 
	MPI_Aint displacements[4];
	MPI_Aint intex, doublex;

	displacements[0] 	= offsetof(simple_c, ll);
	displacements[1] 	= offsetof(simple_c, noise_ll);
	displacements[2] 	= offsetof(simple_c, IDS);
	displacements[3] 	= offsetof(simple_c, ps);
	
	

	MPI_Type_create_struct( 4, blocklens, displacements, old_types, &mystruct );
	MPI_Type_commit( &mystruct );
	
	for (int j = 0; j < root_simple_fits.size(); j++){
		MPI_Send(&root_simple_fits[j], 1, mystruct, 0, j, MPI_COMM_WORLD);
	}




}

void gather_all_segments(){ //for root

}

void send_all_segments(){ //slaves

}
struct bounds{
public:
	double lower_upper[2];
};
struct bounds2{
public:
	int D[3];
};
map<string , vector<vector<double> > > gather_all_bidir_predicitions(vector<segment *> all, 
	vector<segment *> segments , 
	int rank, int nprocs, string out_file_dir, string job_name, int job_ID, params * P, ofstream& log_file){

	map<string , vector<vector<double> > > G;
	map<string , vector<vector<double> > > A;
	//insert data from root
	int N 	= all.size();
	int S;
	int count 	= N/ nprocs;
	vector<bounds2> collections;
	vector<bounds2> final_collections;
	
	if (count==0){
		count 	= 1;
	}
	bounds B;
	bounds2 B2;
	MPI_Datatype mystruct;
	MPI_Datatype mystruct2;
	
	int blocklens[1]={2};
	int blocklens2[1]={3};
	
	MPI_Datatype old_types[1] = {MPI_DOUBLE}; 
	MPI_Datatype old_types2[1] = {MPI_INT}; 
	MPI_Aint displacements[1];
	MPI_Aint displacements2[1];
	
	displacements[0] 	= offsetof(bounds, lower_upper);
	displacements2[0] 	= offsetof(bounds2, D);
	
	MPI_Type_create_struct( 1, blocklens, displacements, old_types, &mystruct );
	MPI_Type_commit( &mystruct );
	
	MPI_Type_create_struct( 1, blocklens2, displacements2, old_types2, &mystruct2 );
	MPI_Type_commit( &mystruct2 );
	if (rank == 0){
		for (int i = 0 ; i < segments.size(); i++){
			for (int j = 0; j < segments[i]->bidirectional_bounds.size();j++){

				vector<double> BBB(2);
				BBB[0]=segments[i]->bidirectional_bounds[j][0], BBB[1]=segments[i]->bidirectional_bounds[j][1];
				G[segments[i]->chrom].push_back(BBB);
				bounds2 B2;
				B2.D[0]=i, B2.D[1]=int(BBB[0]),B2.D[2]=int(BBB[1]);
				collections.push_back(B2);
			}
		}
		//count 	= int(collections.size()) / nprocs;
		for (int j =1; j < nprocs; j++){
		  int start 	= j * count;
		  int stop 	= min(start + count, int(N ));	
	     
		  if (j==nprocs-1){
		    stop 	= N;
		  }
		  if (start >= stop){
		    start 	= stop;
		  }
		 
		 
		  
		  
		  for (int i = 0; i < (stop-start) ; i++){
				MPI_Recv(&S, 1, MPI_INT, j, i, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			     
				
				//MPI_Barrier(MPI_COMM_WORLD);
				G[all[ start+i ]->chrom] 	= vector<vector<double> >(S);
				for (int u = 0; u < S; u++ ){
				 
					MPI_Recv(&B, 1, mystruct, j, u, MPI_COMM_WORLD,MPI_STATUS_IGNORE);			
				       
					//	MPI_Barrier(MPI_COMM_WORLD);
					G[all[ start+i ]->chrom][u] 	= {B.lower_upper[0],B.lower_upper[1]};
					bounds2 B2;
					B2.D[0]=start+i, B2.D[1]=int(B.lower_upper[0]),B2.D[2]=int(B.lower_upper[1]);
					collections.push_back(B2);
				}
			}
		}
		string OUT_ID 	= "Total Number of moment estimator predictions: "+ to_string(int(collections.size()));
		printf("%s\n", OUT_ID.c_str());
		log_file<<OUT_ID<<endl;


	}else  {
	  int IS = segments.size();
	  
		for (int i = 0;  i < segments.size();i++ ){
			int S 	= segments[i]->bidirectional_bounds.size();
			
			MPI_Ssend(&S, 1, MPI_INT, 0,i,MPI_COMM_WORLD );
			//MPI_Barrier(MPI_COMM_WORLD);
			
			for (int u=0; u < segments[i]->bidirectional_bounds.size(); u++){				
				bounds B;
				B.lower_upper[0] 		= segments[i]->bidirectional_bounds[u][0];
				B.lower_upper[1] 		= segments[i]->bidirectional_bounds[u][1];				
				
				MPI_Send(&B, 1, mystruct, 0, u, MPI_COMM_WORLD);
				//MPI_Barrier(MPI_COMM_WORLD);
				
			}
		}
	}
	if (rank==0 and not out_file_dir.empty()){
		write_out_bidirs(G, out_file_dir, job_name, job_ID, P);
	}
	if (rank==0){
		N 	= collections.size();

		int count 	= N/ nprocs;
		if (count==0){
			count 	= 1;
		}
		for (int j = 1; j < nprocs; j++){
			int start 	= j*count;
			int stop 	= min(start + count, int(N ));	
			if (j==nprocs-1){
				stop 	= N;
			}
			if (start >= stop){
		    	start 	= stop;
		    }
		    
			int S 	= stop-start;
			MPI_Send(&S, 1, MPI_INT, j,1, MPI_COMM_WORLD);
			//MPI_Barrier(MPI_COMM_WORLD);
			//okay now get ready! 
			for (int i = 0; i < S;i++){
				MPI_Send(&collections[start+i], 1, mystruct2, j,i,MPI_COMM_WORLD);
				//	MPI_Barrier(MPI_COMM_WORLD);
			}
		}
		for (int i = 0; i < count;i++){
			final_collections.push_back(collections[i]);
		}
		
	}else{
		MPI_Recv(&S, 1, MPI_INT, 0, 1, MPI_COMM_WORLD,MPI_STATUS_IGNORE);	
		//MPI_Barrier(MPI_COMM_WORLD);
		for (int i = 0; i < S; i++){
		  	MPI_Recv(&B2, 1, mystruct2, 0,i,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			//	MPI_Barrier(MPI_COMM_WORLD);
			final_collections.push_back(B2);
		}
	}
	char processor_name[MPI_MAX_PROCESSOR_NAME];
	int namelen;  
	MPI_Get_processor_name(processor_name, &namelen); 
	for (int i = 0 ; i < final_collections.size(); i++){
		vector <double> BB(2);
		BB[0] = final_collections[i].D[1], BB[1]= final_collections[i].D[2];
		string current_chrom(all[final_collections[i].D[0]]->chrom);
		A[current_chrom ].push_back(BB);
	}
	return A;
}


rsimple_c::rsimple_c(){
};
rsimple_c transform(simple_c sc, vector<segment *> segments ){
	rsimple_c rc;
	rc.st_sp[0]=segments[sc.IDS[0]]->start,rc.st_sp[1]=segments[sc.IDS[0]]->stop;
	rc.st_sp[2]=sc.IDS[0], 	rc.st_sp[3]=sc.IDS[1], rc.st_sp[4]=sc.IDS[2];
	rc.st_sp[5]=segments[sc.IDS[0]]->chrom_ID;
	rc.ps[0]=sc.noise_ll,rc.ps[1]=sc.ll;
	for (int i = 0; i < 13; i++){
		rc.ps[i+2]=sc.ps[i];
	}
	return rc;
}

map<string, map<int, vector<rsimple_c> > > gather_all_simple_c_fits(vector<segment *> segments, 
	vector<simple_c> fits, int rank, int nprocs,map<int, string> chromosomes)
{

	//make MPI derived data type
	rsimple_c rc;
	MPI_Datatype mystruct;
	
	int blocklens[2]={6,15};
	MPI_Datatype old_types[2] = {MPI_INT, MPI_DOUBLE}; 
	MPI_Aint displacements[2];
	displacements[0] 	= offsetof(rsimple_c, st_sp);
	displacements[1] 	= offsetof(rsimple_c, ps);
	
	
	MPI_Type_create_struct( 2, blocklens, displacements, old_types, &mystruct );
	MPI_Type_commit( &mystruct );
	
	vector<rsimple_c> rsimple_c_fits;
	for (int i = 0; i < fits.size(); i++){
		rsimple_c_fits.push_back(transform(fits[i], segments ));
	}
		//first want to send the number of bidirs 
	int S;
	MPI_Status stat; 
	if (rank==0){
		for (int j = 1; j < nprocs; j++){
		 
			MPI_Recv(&S, 1, MPI_INT, j, 0, MPI_COMM_WORLD,&stat);
			
			for (int b = 0; b < S; b++){
			  	MPI_Recv(&rc, 1, mystruct, j, b+1, MPI_COMM_WORLD,&stat);					
			       	rsimple_c_fits.push_back(rc);
			}
		}				
	}else{	
        S 	= int(rsimple_c_fits.size());
	
		MPI_Send(&S, 1, MPI_INT, 0,0, MPI_COMM_WORLD);
		
		for (int b = 0; b < S; b++){
		  
		  
		  	MPI_Send(&rsimple_c_fits[b], 1, mystruct, 0, b+1, MPI_COMM_WORLD);			
		  			  		 
		}
	}
	map<string, map<int, vector<rsimple_c> > > G;
	if (rank==0){
		//want to reorient
		int K;
		string ID;
		for (int i =0; i < rsimple_c_fits.size(); i++)
		{
			rc 	= rsimple_c_fits[i];
			string chrom 	= chromosomes[rc.st_sp[5]];
			ID 	= string(chrom) + ":" + to_string(rc.st_sp[0]) + "-" + to_string(rc.st_sp[1]);
			K 	= rc.st_sp[3] ;
			G[ID][K].push_back(rc);
		}
		
	}


	return G;
}


struct seg_and_bidir{
	char chrom[5];
	int st_sp[4];
	double parameters[4];
	seg_and_bidir(segment * current, vector<double> p){
		for (int c = 0; c < 5; c++){
			if (c < current->chrom.size()){
				chrom[c] 	= current->chrom[c];
			}else{
				chrom[c] 	= '\0';
			}
		}
		st_sp[0] 		= current->start,st_sp[1] 	= current->stop;
		st_sp[2] 		= current->ID;
		parameters[0] 	= p[2],parameters[1] 	= p[3];
		parameters[2] 	= p[4],parameters[3] 	= p[5];
		
	}
	seg_and_bidir(){};
};

seg_and_bidir seg_to_seg_and_bidir(segment * s, vector<double> ps, int i ){
	seg_and_bidir sb;
	for (int i = 0; i < 5; i++){
		if (i < s->chrom.size()){
			sb.chrom[i] 	= s->chrom[i];
		}else{
			sb.chrom[i] 	= '\0';
		}
	}
	sb.st_sp[0] 		= i, sb.st_sp[1] = s->start, sb.st_sp[2] = s->stop, sb.st_sp[3]=s->ID;
	sb.parameters[0] 	= ps[2],sb.parameters[1] 	= ps[3];
	sb.parameters[2] 	= ps[4],sb.parameters[3] 	= ps[5];
	return sb;
	
}


struct simple_seg_struct{
	char chrom[6];
	char strand[2];
	int st_sp[4]; //first->start, second->stop
};

vector<segment *> segment_sort(vector<segment *> X){
	bool changed=true;
	while (changed){
		changed=false;
		for (int i = 0; i < X.size()-1; i++  )	{
			if (X[i]->start > X[i+1]->start){ //sort by starting position
				segment * copy 			= X[i];
				X[i] 					= X[i+1];
				X[i+1] 					= copy;
				changed=true;
			}
		}
	}
	return X;
}


map<string, vector<segment *> > send_out_elongation_assignments(vector<segment *> FSI, int rank, int nprocs){
	seg_and_bidir sb;
	MPI_Datatype mystruct;
	
	int blocklens[3]={5,4,4};
	
	MPI_Datatype old_types[3] = {MPI_CHAR, MPI_INT, MPI_DOUBLE}; 
	
	MPI_Aint displacements[3];
	
	displacements[0] = offsetof(seg_and_bidir, chrom),displacements[1] = offsetof(seg_and_bidir, st_sp),displacements[2] = offsetof(seg_and_bidir, parameters);
	
	MPI_Type_create_struct( 3, blocklens, displacements, old_types, &mystruct );
	
	MPI_Type_commit( &mystruct );
	
	map<string, vector<segment *> > 	final_out;
	vector<seg_and_bidir> recieved_sb; 
	if (rank== 0){
		//each segment * in FSI has a variable number of bidir predictions...
		//need to send each and keep track...
	
		//1. first make a vector<> of simple_seg_struct
		vector<seg_and_bidir> send_outs;
		typedef vector<segment *>::iterator FSI_type;
		typedef vector<vector<double> >::iterator fitted_type;
		for (FSI_type f = FSI.begin(); f!=FSI.end(); f++){
			for (fitted_type b = (*f)->fitted_bidirs.begin(); b != (*f)->fitted_bidirs.end(); b++ )	{
				seg_and_bidir sab( (*f), (*b)  );
				send_outs.push_back(sab);
			}
		}

		//now we need to send these out to all the other slave nodes
		int N 		= send_outs.size();
		int counts  = N / nprocs;
		if (counts == 0){
			counts	= 1;
		}
		map<int, vector<int> > assignments; 
		typedef vector<seg_and_bidir>::iterator sab_type;
		int prev_ID 		= -100000000;
		int start 			= 0;
		int stop 			= 0;
		int counter 		= 0;
		int j 				= 0;
		for (sab_type s 	= send_outs.begin(); s!= send_outs.end(); s++ ){
			if ((*s).st_sp[2]!=prev_ID){
				if (counter>=counts and j+1 < nprocs){
					vector<int > bounds(3);
					bounds[0] 	= start, bounds[1] 	= stop, bounds[2]= stop-start;
					assignments[j] 	= bounds;
					start 	= stop;
					counter 	= 0;
					j++;
				}
			}
			prev_ID 		= (*s).st_sp[2];
			counter++;
			stop++;
		}
		vector<int > bounds(3);
		bounds[0] 	= start, bounds[1] 	= stop, bounds[2]= stop-start;
		assignments[j] 	= bounds;
		int S 	= 0;
		for (int j = 1; j < nprocs; j++){
			S  	= 0;
			if (assignments.find(j) != assignments.end() ){
				S 	=  assignments[j][2];
			}
			MPI_Send(&S, 1, MPI_INT, j,1, MPI_COMM_WORLD);
			if (S!=0){
				int u 	= 0;
				for (int i = assignments[j][0]; i < assignments[j][1]; i++ ){
					MPI_Send(&send_outs[i], 1, mystruct, j, u, MPI_COMM_WORLD);	
					u++;
				}
			}
		}
		j 	= 0;
		if (assignments.find(j) != assignments.end()){
			for (int i = assignments[j][0]; i < assignments[j][1]; i++ ){
				recieved_sb.push_back(send_outs[i]);
			}	
		}




	}else{
		int S;
		MPI_Recv(&S, 1, MPI_INT, 0, 1, MPI_COMM_WORLD,MPI_STATUS_IGNORE);	
		for (int i 	= 0; i < S; i++){
			MPI_Recv(&sb, 1, mystruct, 0,i,MPI_COMM_WORLD, MPI_STATUS_IGNORE);			
			recieved_sb.push_back(sb);
		}
	}
	map<int, segment * > almost_there;
	typedef vector<seg_and_bidir>::iterator  rsb_type;
	for (rsb_type r = recieved_sb.begin(); r!= recieved_sb.end(); r++ )	{
		if (almost_there.find( (*r).st_sp[2] ) == almost_there.end()  ){
			almost_there[(*r).st_sp[2]] 	= new segment((*r).chrom, (*r).st_sp[0], (*r).st_sp[1],(*r).st_sp[2] );
		}
		vector<double> fb 	= {(*r).parameters[0],(*r).parameters[1], (*r).parameters[2], (*r).parameters[3] };
		almost_there[(*r).st_sp[2]]->add_fitted_bidir(fb);
	}

	//now finally to final_out
	typedef map<int, segment * >::iterator at_type;
	for (at_type a = almost_there.begin(); a!= almost_there.end(); a++ ){
		final_out[a->second->chrom].push_back(a->second);
	}





	return final_out;

}

map<string, vector<segment *> > send_out_single_fit_assignments(vector<segment *> FSI, int rank, int nprocs ){
	map<string, vector<segment *> > GG;
	int N 		= FSI.size();
	int count 	= N / nprocs;
	if (count == 0){
		count 	= 1;
	}
	simple_seg_struct sss;
	MPI_Datatype mystruct;
	
	int blocklens[3]={6,2,4};
	MPI_Datatype old_types[3] = {MPI_CHAR,MPI_CHAR, MPI_INT}; 
	MPI_Aint displacements[3];
	displacements[0] 	= offsetof(simple_seg_struct, chrom);
	displacements[1] 	= offsetof(simple_seg_struct, strand);
	displacements[2] 	= offsetof(simple_seg_struct, st_sp);
	
	
	MPI_Type_create_struct( 3, blocklens, displacements, old_types, &mystruct );
	MPI_Type_commit( &mystruct );

	vector<simple_seg_struct> runs;
	int start, stop;
	int S;
	if (rank == 0){
		//first send out the number you are going to send
		for (int j =0; j < nprocs; j++){
			start 	= j*count;
			stop 	= (j+1)*count;
			stop 		= min(stop, N);
			if (j == nprocs-1){
				stop 	= N;
			}
			if (start > stop){
				stop 	= start; //in this case there are more nodes than segments...so
			}
			S 		= stop-start;
			if (j > 0){
				MPI_Send(&S, 1, MPI_INT, j,1, MPI_COMM_WORLD);
			}
			//now send out structs
			int u 	= 0;
			for (int i =start; i< stop; i++){
				simple_seg_struct SSS;
				for (int c = 0; c < 6; c++){
					if (c <FSI[i]->chrom.size() ){
						SSS.chrom[c] 	= FSI[i]->chrom[c];
					}else{
						SSS.chrom[c] 	= '\0';
					}
				}
				SSS.chrom[5]	= '\0';
				SSS.strand[0] 	= FSI[i]->strand[0];
				SSS.strand[1] 	= '\0';
				SSS.st_sp[0] 	= FSI[i]->start;
				SSS.st_sp[1] 	= FSI[i]->stop;
				SSS.st_sp[2] 	= FSI[i]->ID;
				SSS.st_sp[3] 	= FSI[i]->counts;
				if (j >0){
					MPI_Send(&SSS, 2, mystruct, j, u, MPI_COMM_WORLD  );
				}else{
					runs.push_back(SSS);
				}
				u++;
			}
		}

	}else{
		MPI_Recv(&S, 1, MPI_INT, 0, 1, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		for (int u = 0; u < S; u++){
			MPI_Recv(&sss, 2, mystruct, 0,u,MPI_COMM_WORLD, MPI_STATUS_IGNORE);	
			runs.push_back(sss);
		}	
	}
	//now convert to GG type
	for (int i = 0; i < runs.size(); i++){
		segment * ns 	= new segment(runs[i].chrom, runs[i].st_sp[0], runs[i].st_sp[1], runs[i].st_sp[2], runs[i].strand);
		ns->counts 		= runs[i].st_sp[3];
		GG[ns->chrom].push_back(ns) ;
	}

	typedef map<string, vector<segment *> >::iterator it_type;


	return GG;
}


vector<single_simple_c> gather_all_simple_c(vector<single_simple_c> fits , int rank, int nprocs ){
	single_simple_c sc;
	MPI_Datatype mystruct;
	
	int blocklens[3]={6,3, 10};
	MPI_Datatype old_types[3] = {MPI_CHAR, MPI_INT, MPI_DOUBLE}; 
	MPI_Aint displacements[3];
	displacements[0] 	= offsetof(single_simple_c, chrom);
	displacements[1] 	= offsetof(single_simple_c, st_sp);
	displacements[2] 	= offsetof(single_simple_c, ps);
	
	
	MPI_Type_create_struct( 3, blocklens, displacements, old_types, &mystruct );
	MPI_Type_commit( &mystruct );

	int S;
	vector<single_simple_c> recieved;
	if (rank==0){
		for (int j = 0; j < nprocs; j++){
			if (j > 0){
				MPI_Recv(&S, 1, MPI_INT, j, 1, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				for (int b = 0; b < S; b++){
					MPI_Recv(&sc, 1, mystruct, j, b, MPI_COMM_WORLD,MPI_STATUS_IGNORE);	
					recieved.push_back(sc);
				}
			}else{
				for (int u = 0; u < fits.size(); u++){
					recieved.push_back(fits[u]);
				}
			}
			
		}

	}else{
		//first send the number fitted components
		S 	= fits.size();
		MPI_Send(&S, 1, MPI_INT, 0,1, MPI_COMM_WORLD);
		for (int u = 0; u < S; u++)	{
			MPI_Send(&fits[u], 1, mystruct, 0, u, MPI_COMM_WORLD);
		}
	}
	return recieved;

}

int get_job_ID(string path, string job_ID, int rank, int nprocs){
	//string OUT = dir+"EMGU-" + to_string(job_ID) +"_" + DT+ ".log";
	//tmp_EMGU-0_2.log
	int current 			= 1;
	string current_file 	= "";
	vector<string> line_array;
	int jN 	= job_ID.size();
	bool FOUND = false;
	if (rank==0){
		DIR * dirFile = opendir( path.c_str() );
		if ( dirFile ){
			struct dirent* hFile;
			errno = 0;
			while (( hFile = readdir( dirFile )) != NULL ){
				current_file 	= hFile->d_name;
				if (current_file.substr(0,jN) == job_ID){
					if (current_file.substr(jN, 1) =="-"){
						int t 		= 1;
						while (current_file.substr(jN+t,1) != "_" and t < current_file.size() ){
							t++;
						}
						if (t < current_file.size()){
							current 	= max(current, stoi( current_file.substr(jN+1,t) ) );	
						}else{
							printf("WHOOOO, MPI_comm getting job_ID\n");

						}
						FOUND 		= true;
					}
				}
				
				if (current_file.substr(0, jN + 4 ) =="tmp_" +job_ID){

					if (current_file.substr(jN+4,1) =="-"){
						int t 		= 1;
						while (current_file.substr(jN+4+t,1)   != "_" and t < current_file.size() ){
							t++;
						}
						if (t < current_file.size()){
							current 	= max(current, stoi( current_file.substr(jN+5,t) ) );
						}else{
							printf("WHOOOO, MPI_comm getting job_ID\n");
						}
						FOUND 		= true;
					}
					
				}

			}	
			closedir( dirFile );
		}else{
			cout<<"couldn't open temp log directory"<<endl;
		}
		if (FOUND){
			current++;
		}
		for (int j =1; j < nprocs; j++ ){
			MPI_Send(&current, 1, MPI_INT, j,1, MPI_COMM_WORLD);		
		}


	}else{
		MPI_Recv(&current, 1, MPI_INT, 0, 1, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		
	}
	return current;
}

	
map<int, map<int, vector<simple_c_free_mode>  > > gather_all_simple_c_free_mode(vector<map<int, vector<simple_c_free_mode> >> FITS, 
	int rank, int nprocs, ofstream& FHW){
	

	simple_c_free_mode sc_fm;
	MPI_Datatype mystruct;
	
	int blocklens[4]={3,5, 6, 12};
	MPI_Datatype old_types[4] = {MPI_DOUBLE, MPI_INT, MPI_CHAR, MPI_DOUBLE}; 
	MPI_Aint displacements[4];
	displacements[0] 	= offsetof(simple_c_free_mode, SS);
	displacements[1] 	= offsetof(simple_c_free_mode, ID);
	displacements[2] 	= offsetof(simple_c_free_mode, chrom);
	displacements[3] 	= offsetof(simple_c_free_mode, ps);
	
	
	MPI_Type_create_struct( 4, blocklens, displacements, old_types, &mystruct );
	MPI_Type_commit( &mystruct );


	vector<simple_c_free_mode> recieved;

	typedef map<int, vector<simple_c_free_mode> >::iterator model_it;
	if (rank==0){
		//first want to recieve the number of segments each child processes ran
		for (int j = 0; j < nprocs; j++){
			int S =0;
			if (j > 0){
				FHW<<"(MPI_comm; root), waiting on " + to_string(j) +"\n";
				FHW.flush();
				MPI_Recv(&S, 1, MPI_INT, j, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				FHW<<"(MPI_comm; root), recieved from " + to_string(j) + ", " + to_string(S) + " fits\n";
				FHW.flush();
				FHW<<"(MPI_comm; root), waiting on " + to_string(j) + ", to receive all " + to_string(S) + " fits\n";
				FHW.flush();
				for (int s = 0; s < S; s++){
					FHW<<to_string(s+1)<<":";
					FHW.flush();
					MPI_Recv(&sc_fm, 1, mystruct,j,s+1,  MPI_COMM_WORLD,MPI_STATUS_IGNORE );
					FHW<<to_string(s+1)<<",";
					FHW.flush();
					recieved.push_back(sc_fm);
				}
				FHW<<"(MPI_comm; root), recieved from " + to_string(j) + ", all " + to_string(S) + " fits\n";
				
			}else{
				for (int s = 0; s < FITS.size(); s++){
					for (model_it k = FITS[s].begin(); k!=FITS[s].end(); k++){
						for (int sc= 0; sc < k->second.size();sc++ ){
							recieved.push_back( k->second[sc] );
						}
					}
				}					
			}
		}

	}else{

		int S 	= 0;
		for (int s = 0; s < FITS.size(); s++){
			for (model_it k = FITS[s].begin(); k!=FITS[s].end(); k++){
				for (int sc= 0; sc < k->second.size();sc++ ){
					S++;
				}
			}
		}
		FHW<<"(MPI_comm; slave process: " + to_string(rank)+ " about to send " + to_string(S) + " fits\n";
		FHW.flush();
		MPI_Send(&S, 1, MPI_INT, 0,0, MPI_COMM_WORLD);
		FHW<<"(MPI_comm; slave process: " + to_string(rank)+ " sent " + to_string(S) + " fits\n";
		FHW.flush();
		S=0;
		FHW<<"(MPI_comm; slave process: " + to_string(rank)+ " sending all " + to_string(S) + " fits\n";
		FHW.flush();
		for (int s = 0; s < FITS.size(); s++){
			for (model_it k = FITS[s].begin(); k!=FITS[s].end(); k++){
				for (int sc= 0; sc < k->second.size();sc++ ){

					MPI_Send(&k->second[sc], 1, mystruct, 0, S+1, MPI_COMM_WORLD);	
					FHW<<to_string(S+1)<<",";
					S++;
				}
			}
		}
		FHW<<"(MPI_comm; slave process: " + to_string(rank)+ " sent all " + to_string(S) + " fits\n";
		FHW.flush();
		
	}
	//printf("Rank: %d,%d\n",rank, recieved.size());
	map<int, map<int, vector<simple_c_free_mode>  > > G;
	typedef vector<simple_c_free_mode>::iterator it_type_fm;

	for (it_type_fm sc = recieved.begin(); sc!=recieved.end(); sc++){
		G[(*sc).ID[0]][(*sc).ID[3]].push_back(*sc);
	}

	return G;



}
struct st_sp_bidir{
	int st_sp[2];
};

boostrap_struct::boostrap_struct(){}
boostrap_struct::boostrap_struct(vector<double> old, 
	vector<double> variances,int start, int stop, int chrom_ID ){
		IDS[0] 	= start, IDS[1] = stop, IDS[2] = chrom_ID;
		int j 	= 0;
		for (int i = 0; i < old.size(); i++){
			parameters[j] 	= old[i];
			j++;
		}
		for (int i = 0; i < variances.size(); i++){
			parameters[j] 	= variances[i];
			j++;
		}
}
string boostrap_struct::print_out(map<int, string> G, params * P){
	double scale 	= stod(P->p6["-ns"]);
	string line ="",old="", var="";
	if (G.find(IDS[2])==G.end()  ){
		printf("Chromosome not in ID_to_chrom...bootsrap_struct...\n");
		return line;
	}
	int center 	= parameters[0] ;
	int std 	= (parameters[1] /2.)  + ( parameters[2]);
	line 		= G[IDS[2]] + "\t" + to_string(center-std) + "\t" + to_string(center+std)+ "\t";
	for (int i = 0 ; i < 7; i++){
		old+=to_string(parameters[i])+"_";
	}
	old+=to_string(parameters[7]) + "\t";

	for (int i=8; i < 13;i++){
		var+=to_string(parameters[i])+"_";
	}
	var+=to_string(parameters[13]);
	return line + old+ var+ "\n";

}

vector<int> send_out_merged_start_stops(vector<vector<int>> start_stops, int rank, int nprocs){
	st_sp_bidir ssb;
	MPI_Datatype mystruct;
	
	int blocklens[1]={2};
	MPI_Datatype old_types[1] = { MPI_INT }; 
	MPI_Aint displacements[1];
	displacements[0] 	= offsetof(st_sp_bidir, st_sp);
	
	
	MPI_Type_create_struct( 1, blocklens, displacements, old_types, &mystruct );
	MPI_Type_commit( &mystruct );

	
	vector<int> st_sp(2) ;
	if (rank==0){
		st_sp = start_stops[rank];
		for (int j = 1 ; j < nprocs; j++){
			st_sp_bidir ssb2;

			ssb2.st_sp[0] 	= start_stops[j][0],ssb2.st_sp[1] 	= start_stops[j][1];
			MPI_Send(&ssb2, 1, mystruct, j, 0, MPI_COMM_WORLD);

		}
	}else{
		MPI_Recv(&ssb, 1, mystruct, 0, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);	
		st_sp[0] 	= ssb.st_sp[0], st_sp[1]=ssb.st_sp[1];
	}
	return st_sp;
}

vector<boostrap_struct> collect_bootstrap(vector<segment *> segments, int rank, int nprocs,
	map<string, int> chrom_to_ID){
	//need to send 
	//chrom,start, stop -> ints
	//MLE paramters from before
	//variances of these parameters
	//old LL, and old N
	boostrap_struct bs;
	MPI_Datatype mystruct;
	
	int blocklens[2]={3,14};
	MPI_Datatype old_types[4] = { MPI_INT, MPI_DOUBLE}; 
	MPI_Aint displacements[2];
	displacements[0] 	= offsetof(boostrap_struct, IDS);
	displacements[1] 	= offsetof(boostrap_struct, parameters);
	
	
	MPI_Type_create_struct( 2, blocklens, displacements, old_types, &mystruct );
	MPI_Type_commit( &mystruct );


	int S ;
	vector<boostrap_struct> BS;
	for (int i = 0 ; i < segments.size(); i++){
		if (segments[i]->parameters.size() == segments[i]->variances.size()) {
			for (int j = 0; j < segments[i]->parameters.size(); j++){
				BS.push_back(boostrap_struct(segments[i]->parameters[j], segments[i]->variances[j], 
					segments[i]->start, segments[i]->stop, chrom_to_ID[segments[i]->chrom] )  );
			}
		}else{
			printf("Variances != Old Parameters Size....MPI_comm\n");
		}
	}
	if (rank==0){
		for (int j =1 ; j < nprocs; j++){
			MPI_Recv(&S, 1, MPI_INT, j, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			for (int s = 0; s < S; s++)	{
				MPI_Recv(&bs, 1, mystruct,j,s+1,  MPI_COMM_WORLD,MPI_STATUS_IGNORE );
				BS.push_back(bs);
			}
		}
	}else{
		S 	= BS.size();
		MPI_Send(&S, 1, MPI_INT, 0,0, MPI_COMM_WORLD);
		for (int s= 0; s < S;s++ ){
			MPI_Send(&BS[s], 1, mystruct, 0, s+1, MPI_COMM_WORLD);	
		}
	}
	return BS;
}









