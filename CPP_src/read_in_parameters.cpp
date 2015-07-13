#include <iostream>
#include <fstream>
#include <unistd.h>
#include <string>
#include <vector>
#include <algorithm>
#include "read_in_parameters.h"
#include <stdio.h>
#include <ctype.h>
using namespace std;
params::params(){
	p["-i"] 		= "";
	p["-o"] 		= "";
	p["-br"] 		= "300";
	p["-ns"] 		= "100";
	p["-maxK"] 		= "3";
	p["-rounds"] 	= "10";
	p["-ct"] 		= "0.0001";
	p["-max_noise"] = "0.05";
	p["-move"] 		= "5";
	p["-np"] 		= "4";
	p["-chr"] 		= "all";
	p["-v"] 		= "0";
	p["-mi"] 		= "300";
	p["-pp_r"] 		= "1";
	p["-pp_k"] 		= "0";
	p["-r_mu"] 		= "0";
	N 				= 0;

}
void params::help(){
	cout<<"----------------------------------------------------------------"<<endl;
	cout<<"               Description of EMGU Parameters"<<endl;
	cout<<"-i  : "<<"/path/to/annotation_file"<<endl;
	cout<<"-o  : "<<"path/to/output_file"<<endl;
	cout<<"-np : "<<"number of CPU cores (default: 1)"<<endl;
	cout<<"-chr: "<<"specific chromsome (default: all)"<<endl;
	cout<<"-r  : "<<"base pair binning resolution (default: 1)"<<endl;
	cout<<"Questions/Bugs? joseph[dot]azofeifa[at]colorado[dot]edu"<<endl;
	cout<<"These parameters can also be set in seperate config file"<<endl;
	cout<<"----------------------------------------------------------------"<<endl;
	
}
void params::display(){
	cout<<"----------------------------------------------------------------"<<endl;
	cout<<"              User Provided EMGU Parameters                     "<<endl;
	cout<<"-i        : "<<p["-i"]<<endl;
	cout<<"-o        : "<<p["-o"]<<endl;
	cout<<"-chr      : "<<p["-chr"]<<endl;
	cout<<"-br       : "<<p["-br"]<<endl;
	cout<<"-ns       : "<<p["-ns"]<<endl;
	cout<<"-maxK     : "<<p["-maxK"]<<endl;
	cout<<"-rounds   : "<<p["-rounds"]<<endl;
	cout<<"-ct       : "<<p["-ct"]<<endl;
	cout<<"-mi       : "<<p["-mi"]<<endl;
	cout<<"-max_noise: "<<p["-max_noise"]<<endl;
	cout<<"-move     : "<<p["-move"]<<endl;
	cout<<"-np       : "<<p["-np"]<<endl;
	cout<<"-pp_r     : "<<p["-pp_r"]<<endl;
	cout<<"-pp_k     : "<<p["-pp_k"]<<endl;
	cout<<"-r_mu     : "<<p["-r_mu"]<<endl;
	
	cout<<"----------------------------------------------------------------"<<endl;
}

bool checkIfFileAndConfigFile(string FILE){
	ifstream FH(FILE);
	
	if (FH){
		string line;
		getline(FH, line);
		if ("#" == line.substr(0,1)){
			FH.close();
			return true;
		}
		return false;
	}
	return false;
}

void read_in_config_file(string FILE, params * P){
	ifstream FH(FILE);
	if (FH){
		string line, ID, val;
		string::iterator it;
		string empty_space 	= " ";
		string dash  		= "-";
		string equal 		= "=";
		string pound 		= "#";
		
		const char *es 	=empty_space.c_str();
		const char *d 	=dash.c_str();
		const char *eq 	=equal.c_str();
		const char *h 	=pound.c_str();
		
		bool add_ID 	= true;
		bool add_param 	= false;
		while (getline(FH, line)){
			if ("#" != line.substr(0,1)){
				ID 	= "";
				val = "";
				add_ID=true;
				add_param=false;
				for (int i = 0; i < line.size(); i++){
					if (line[i]==*d){
						ID+="-";
					}else if (not ID.empty() and not isspace(line[i]) and add_ID ){
						ID+=line[i];
					}else if (not ID.empty() and isspace(line[i])){
						add_ID=false;
					}
					else if (not ID.empty() and line[i]==*eq ){
						add_param=true;
					}else if (add_param and not isspace(line[i]) and line[i]!=*h){
						val+=line[i];
					}else if (add_param and not isspace(line[i])){
						add_param=false;
					}
				}
				if (P->p.find(ID) !=P->p.end()){
					P->p[ID] 	= val;
				}				
			}	
		}
	}else{
		printf("couldn't open config file...\n");
	}
}


void fillInOptions(char* argv[],params * P){
	string F 		= "";
	char * COM 		= "-";
	bool begin 		= true;
	bool GO_FORIT 	= false;
	while (*argv){
		if (begin){
			F 				= string(*argv); 
			bool isConfig 	= checkIfFileAndConfigFile(F);
			if (isConfig){
				//read this in
				read_in_config_file(F, P);
				begin 		= false;

			}else{
				GO_FORIT 	= true;
			}

		}
		if ((*argv)[0] == COM[0] and GO_FORIT){
			F 			= string(*argv); 
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
		GO_FORIT=true;
	}

}




params * readInParameters( char* argv[]){	
	string userModParameter = "";
	params 	* P = new params;
	fillInOptions(argv, P);
	return P;
}


