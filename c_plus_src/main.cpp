//============================================================================
// Name        : main.cpp
// Author      : Joey Azofeifa
// Version     : 1.0
// Description : Main file for running Mixture Model Package
//============================================================================
#include <iostream>
#include "read_in_parameters.h"
#include "load.h"
#include "write.h"
#include <string>
#include <map>
#include <vector>
#include <iostream>
#include "run_fit.h"

using namespace std;

int main(int argc, char* argv[]) {
    params * P 				= new params();
    P               		= readInParameters(argv);
    //=================================================================
    // so here are a list of the parameters that user can specify
    // ("-i") genome coverage file 
    // ("-j") annotation file
    // ("-o") out file
    // ("-r") bin size resolution, this is the width of each bin
    // ("-np") number of processors
    // ("-chr") specific chromosome to run on
    //=================================================================
    if (stoi(P->p["-h"]) or P->N==0){
    	P->help();
    	return 1;
    }



	string bgFile_pos 		= P->p["-i"];
	string bgFile_neg		= P->p["-j"];
	string annotFile 		= P->p["-k"];

	string outFile 			= P->p["-o"];
	int resolution 			= stoi(P->p["-r"]);
	int num_proc 			= stoi(P->p["-np"]);
	int rt 					= stoi(P->p["-rt"]);
	double penality 		= stod(P->p["-pb"]);
	int maxK 				= stoi(P->p["-mk"]);
	bool test 				= stoi(P->p["-t"]);
	bool verbose 			= stoi(P->p["-v"]);
	bool standardize_bool 	= 1;
	bool DMGU 				= stoi(P->p["-dgu"]);
	bool EMGU 				= stoi(P->p["-egu"]);
	int pad 				= stoi(P->p["-p"]);
	string strand 			= P->p["-s"];
	int S;
	if (strand=="all"){
		S=2;
	}else if(strand=="forward"){
		S=0;
	}else{
		S=1;
	}
	if (not DMGU and not EMGU){
		cout<<"mixture modelling regime not specified..."<<endl;
		cout<<"exiting..."<<endl;
		return 0;
	}

	P->display();
	annotations * A = new annotations;
	if (verbose){
		cout<<"reading annotations: ";
		cout<<flush;
	}
	read_annotation(annotFile, resolution, A, pad);//read in annotation file
	if (verbose){
		cout<<flush;
		cout<<"done"<<endl;
	}
	if (verbose){
		cout<<"reading forward strand bedgraph: ";
		cout<<flush;
	}
	read_bedgraph(bgFile_pos, A, 1);//read in positive strand file
	if (verbose){
		cout<<flush;
		cout<<"done"<<endl;
	}
	if (verbose){
		cout<<"reading reverse strand bedgraph: ";
		cout<<flush;
	}
	read_bedgraph(bgFile_neg, A, 2);//read in negative strand file

	if (verbose){
		cout<<flush;
		cout<<"done"<<endl;
	}
	if (verbose){
		cout<<"center and normalizing read data for numerical stability: ";
		cout<<flush;
	}
	
	if (standardize_bool){
		standardize(A);
	}
	if (verbose){
		cout<<flush;
		cout<<"done"<<endl;
	}
	
	
	fit(A, DMGU, EMGU, verbose, S);

	
	write_out_file(outFile, A, test);

	delete(A);
	return 1;
	
}
