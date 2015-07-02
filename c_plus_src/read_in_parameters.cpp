#include <iostream>
#include <fstream>
#include <unistd.h>
#include <string>
#include <vector>
#include <algorithm>
#include "read_in_parameters.h"
using namespace std;
params::params(){
	p["-i"] 	= "";
	p["-j"] 	= "";
	p["-k"] 	= "";
	p["-o"] 	= "";
	p["-r"] 	= "1";
	p["-chr"] 	= "all";
	p["-np"] 	= "1";
	p["-rt"] 	= "20";
	p["-pb"] 	= "1";
	p["-mk"] 	= "5";
	p["-t"] 	= "0";
	p["-s"]		= "all";
	p["-v"] 	= "0";
	p["-dgu"] 	= "0";
	p["-egu"] 	= "0";
	p["-p"] 	= "0";
	p["-h"] 	= "0";
	N 			= 0;

}
void params::help(){
	cout<<"----------------------------------------------------------------"<<endl;
	cout<<"               Description of EMGU Parameters"<<endl;
	cout<<"-i  : "<<"/path/to/forward_strand_bedgraph_file.bedgraph"<<endl;
	cout<<"-j  : "<<"/path/to/reverse_strand_bedgraph_file.bedgraph"<<endl;
	cout<<"-k  : "<<"/path/to/annotation_file"<<endl;
	cout<<"-o  : "<<"path/to/output_file"<<endl;
	cout<<"-np : "<<"number of CPU cores (default: 1)"<<endl;
	cout<<"-chr: "<<"specific chromsome (default: ALL)"<<endl;
	cout<<"-r  : "<<"base pair binning resolution (default: 1)"<<endl;
	cout<<"Questions/Bugs? joseph[dot]azofeifa[at]colorado[dot]edu"<<endl;
	cout<<"----------------------------------------------------------------"<<endl;
	
}
void params::display(){
	if (p["-v"]=="1"){
		cout<<"----------------------------------------------------"<<endl;
		cout<<"                 Running EMGU"<<endl;
		cout<<"Questions/Bugs? joseph[dot]azofeifa[at]colorado[dot]edu"<<endl;
		cout<<"Forward Strand Input File: "<< p["-i"]<<endl;
		cout<<"Reverse Strand Input File: "<< p["-j"]<<endl;
		cout<<"Annotation File          : "<< p["-k"]<<endl;
		cout<<"Output File              : "<< p["-o"]<<endl;
		cout<<"Processors               : "<< p["-np"]<<endl;
		cout<<"Specific Chromosome      : "<< p["-chr"]<<endl;
		cout<<"bp Resolution            : "<< p["-r"]<<endl;
		cout<<"----------------------------------------------------"<<endl;
	}
	if (p["-v"]=="1"){
		if (p["-DGU"]=="1"){
			cout<<"Running Double Geometric - Uniform Model"<<endl;
			cout<<"Warning: Can only support single strand/mixture fit"<<endl;
		}else if (p["-EGU"]=="1"){
			cout<<"Running Exponentially Modified Gaussian - Uniform Model"<<endl;
		}
	}
}

void fillInOptions(char* argv[],params * P){
	string F 		= "";
	char * COM 		= "-";
	while (*argv){
		if ((*argv)[0] == COM[0]){
			F 			= string(*argv); 
			transform(F.begin(), F.end(), F.begin(), ::tolower);
			if (P->p.find(F) !=P->p.end()){
				P->p[F] 	= "1";
				P->N+=1;
			}else{
				F 			= "";
			}
		}
		else if (not F.empty()) {
			if (P->p.find(F) !=P->p.end()){

				P->p[F]=string(*argv);
			}
		}
		argv++;
	}

}




params * readInParameters( char* argv[]){	
	string userModParameter = "";
	params 	* P = new params;
	fillInOptions(argv, P);
	return P;
}


