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
	p["-v"] 		= "1";
	p["-i"] 		= "";
	p["-o"] 		= "";
	p["-br"] 		= "300";
	p["-ns"] 		= "100";
	p["-minK"] 		= "1";
	p["-maxK"] 		= "3";
	p["-rounds"] 	= "10";
	p["-ct"] 		= "0.0001";
	p["-max_noise"] = "0.05";
	p["-np"] 		= "4";
	p["-chr"] 		= "all";
	p["-v"] 		= "0";
	p["-mi"] 		= "300";
	p["-r_mu"] 		= "0";
	p["-ALPHA_0"] 	= "1";
	p["-BETA_0"] 	= "1";
	p["-ALPHA_1"] 	= "1";
	p["-BETA_1"] 	= "1";
	p["-ALPHA_2"] 	= "1";
	p["-ALPHA_3"] 	= "1";
	p["-template"] 	= "1";
	p["-density"] 	= "1000";
	p["-window"] 	= "2000";
	p["-pad"] 		= "5000";
 	

	p2["-v"] 		= "1";
	p2["-i"] 		= "";
	p2["-j"] 		= "";
	p2["-k"] 		= "";
	p2["-o"] 		= "";
	p2["-pad"] 		= "100";

	p3["-v"] 		= "1";
	p3["-i"] 		= "";
	p3["-o"] 		= "";
	p3["-penality"] = "1";
	p3["-to_igv"] 	= "1";
	p3["-to_EMG"] 	= "1";


	p4["-v"] 				= "1";
	p4["-i"] 				= "";
	p4["-j"] 				= "";
	p4["-k"] 				= "";
	p4["-o"] 				= "1";
	p4["-ns"] 				= "100";
	p4["-br"] 				= "50";
	p4["-density"] 			= "1000";
	p4["-window"] 			= "1500";
	p4["-chr"] 				= "all";
	p4["-opt_res"] 			= "5";
	p4["-np"] 				= "4";



	
	N 				= 0;
	module 			= "";
	EXIT 			= 0;

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
	if (module=="MODEL"){
		if (p["-template"]!="1"){
			cout<<"----------------------------------------------------------------"<<endl;
			cout<<"              User Provided EMGU Parameters                     "<<endl;
			cout<<"                 (RUNNING MIXTURE MODEL)                        "<<endl;
			cout<<"-i        : "<<p["-i"]<<endl;
			cout<<"-o        : "<<p["-o"]<<endl;
			cout<<"-chr      : "<<p["-chr"]<<endl;
			cout<<"-br       : "<<p["-br"]<<endl;
			cout<<"-ns       : "<<p["-ns"]<<endl;
			cout<<"-minK     : "<<p["-minK"]<<endl;
			cout<<"-maxK     : "<<p["-maxK"]<<endl;
			cout<<"-rounds   : "<<p["-rounds"]<<endl;
			cout<<"-ct       : "<<p["-ct"]<<endl;
			cout<<"-mi       : "<<p["-mi"]<<endl;
			cout<<"-max_noise: "<<p["-max_noise"]<<endl;
			cout<<"-np       : "<<p["-np"]<<endl;
			cout<<"-r_mu     : "<<p["-r_mu"]<<endl;
			
			cout<<"----------------------------------------------------------------"<<endl;
			cout<<"Questions/Bugs? joseph[dot]azofeifa[at]colorado[dot]edu"<<endl;
		}
		else{
			cout<<"----------------------------------------------------------------"<<endl;
			cout<<"              User Provided EMGU Parameters                     "<<endl;
			cout<<"                 (RUNNING MIXTURE MODEL)                        "<<endl;
			cout<<"            ...coupled to template matching...                  "<<endl;	
			cout<<"-i        : "<<p["-i"]<<endl;
			cout<<"-o        : "<<p["-o"]<<endl;
			cout<<"-chr      : "<<p["-chr"]<<endl;
			cout<<"-br       : "<<p["-br"]<<endl;
			cout<<"-ns       : "<<p["-ns"]<<endl;
			cout<<"-rounds   : "<<p["-rounds"]<<endl;
			cout<<"-ct       : "<<p["-ct"]<<endl;
			cout<<"-mi       : "<<p["-mi"]<<endl;
			cout<<"-np       : "<<p["-np"]<<endl;
			cout<<"-density  : "<<p["-density"]<<endl;
			cout<<"-window   : "<<p["-window"]<<endl;
			cout<<"-pad      : "<<p["-pad"]<<endl;
			cout<<"----------------------------------------------------------------"<<endl;
			cout<<"Questions/Bugs? joseph[dot]azofeifa[at]colorado[dot]edu"<<endl;
		}
	}else if (module=="FORMAT"){
		cout<<"----------------------------------------------------------------"<<endl;
		cout<<"              User Provided EMGU Parameters                     "<<endl;
		cout<<"                    (FORMATING DATA)                        "<<endl;
		cout<<"-i        : "<<p2["-i"]<<endl;
		cout<<"-j        : "<<p2["-j"]<<endl;
		cout<<"-k        : "<<p2["-k"]<<endl;
		cout<<"-o        : "<<p2["-o"]<<endl;
		cout<<"-pad      : "<<p2["-pad"]<<endl;
		cout<<"----------------------------------------------------------------"<<endl;
		cout<<"Questions/Bugs? joseph[dot]azofeifa[at]colorado[dot]edu"<<endl;
	}else if (module=="SELECTION"){
		cout<<"----------------------------------------------------------------"<<endl;
		cout<<"              User Provided EMGU Parameters                     "<<endl;
		cout<<"                    (MODEL SELECTION)                        "<<endl;
		cout<<"-i           : "<<p3["-i"]<<endl;
		cout<<"-o           : "<<p3["-i"]<<endl;
		cout<<"-penality    : "<<p3["-penality"]<<endl;
		cout<<"-to_igv      : "<<p3["-to_igv"]<<endl;
		cout<<"-to_EMG      : "<<p3["-to_EMG"]<<endl;
		cout<<"----------------------------------------------------------------"<<endl;
		cout<<"Questions/Bugs? joseph[dot]azofeifa[at]colorado[dot]edu"<<endl;
	}else if (module=="BIDIR")	{
		if (p4["-k"].empty()){
			cout<<"----------------------------------------------------------------"<<endl;
			cout<<"              User Provided EMGU Parameters                     "<<endl;
			cout<<"                 (BIDIRECTIONAL DETECTOR)                       "<<endl;
			cout<<"-i           : "<<p4["-i"]<<endl;
			cout<<"-j           : "<<p4["-j"]<<endl;
			cout<<"-k           : "<<p4["-k"]<<endl;
			cout<<"-o           : "<<p4["-o"]<<endl;
			cout<<"-window      : "<<p4["-window"]<<endl;
			cout<<"-density     : "<<p4["-density"]<<endl;
			cout<<"-ns          : "<<p4["-ns"]<<endl;
			cout<<"-br          : "<<p4["-br"]<<endl;
			cout<<"-np          : "<<p4["-np"]<<endl;
			cout<<"----------------------------------------------------------------"<<endl;
			cout<<"Questions/Bugs? joseph[dot]azofeifa[at]colorado[dot]edu"<<endl;
		}else{
			cout<<"----------------------------------------------------------------"<<endl;
			cout<<"              User Provided EMGU Parameters                     "<<endl;
			cout<<"                 (BIDIRECTIONAL DETECTOR)                       "<<endl;
			cout<<"                ...Running Optimization...                      "<<endl;
			cout<<"-i           : "<<p4["-i"]<<endl;
			cout<<"-j           : "<<p4["-j"]<<endl;
			cout<<"-k           : "<<p4["-k"]<<endl;
			cout<<"-o           : "<<p4["-o"]<<endl;
			cout<<"-opt_res     : "<<p4["-opt_res"]<<endl;
			cout<<"-chr         : "<<p4["-chr"]<<endl;
			cout<<"-np          : "<<p4["-np"]<<endl;
			cout<<"-ns          : "<<p4["-ns"]<<endl;
			cout<<"-br          : "<<p4["-br"]<<endl;
			cout<<"----------------------------------------------------------------"<<endl;
			cout<<"Questions/Bugs? joseph[dot]azofeifa[at]colorado[dot]edu"<<endl;
			
		}
	}
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

string checkModule(string FILE){
	ifstream FH(FILE);
	
	if (FH){
		string line, ID;
		ID="";
		string pound 		= "#";
		const char *h 	=pound.c_str();
		while (getline(FH, line)){
			if ("~"==line.substr(0,1)){
				for (int i = 1; i < line.size(); i++){
					if (not isspace(line[i]) and line[i]!= *h){
						ID+=line[i];
					}else{
						break;
					}
				}
				FH.close();
				return ID;
			}
		}
	}
	return "";
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
		vector<string> not_valid;
		while (getline(FH, line)){
			if ("#" != line.substr(0,1)){
				ID 	= "";
				val = "";
				add_ID=true;
				add_param=false;
				for (int i = 0; i < line.size(); i++){
					if (line[i]==*d and not add_param and ID.empty()){
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
				if (P->module=="MODEL"){
					if (P->p.find(ID) !=P->p.end()){
						P->p[ID] 	= val;
					}else if(not ID.empty()) {
						not_valid.push_back(ID);
					}				
				}else if (P->module=="FORMAT"){
					if (P->p2.find(ID) !=P->p2.end()){
							P->p2[ID] 	= val;
					}else if (not ID.empty() ){
						not_valid.push_back(ID);
					}
				}else if (P->module=="SELECTION")	{
					if (P->p3.find(ID) !=P->p3.end()){
							P->p3[ID] 	= val;
					}else if (not ID.empty() ){
						not_valid.push_back(ID);
					}
				}else if (P->module=="BIDIR")	{
					if (P->p4.find(ID) !=P->p4.end()){
							P->p4[ID] 	= val;
					}else if (not ID.empty() ){
						not_valid.push_back(ID);
					}
				}
			}	
		}
		if (not_valid.size() >0){
			cout<<"These parameters were found but are not valid identifiers: ";
			for (int i =0; i < not_valid.size(); i++){
				cout<<not_valid[i]<<",";
			}
			P->EXIT=1;
			cout<<endl;
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
				string module 	= checkModule(F);
				P->module 		= module;
				if (not P->module.empty() ){
					read_in_config_file(F, P);
				}
				begin 		= false;
			}else{
				GO_FORIT 	= true;
			}

		}
		if ((*argv)[0] == COM[0] and GO_FORIT){
			F 			= string(*argv); 
			if (P->module=="MODEL"){
				if (P->p.find(F) !=P->p.end()){
					P->p[F] 	= "1";
					P->N+=1;
				}else{
					F 			= "";
				}
			}else if(P->module=="FORMAT"){
				if (P->p2.find(F) !=P->p2.end()){
					P->p2[F] 	= "1";
					P->N+=1;
				}else{
					F 			= "";
				}
			}else if (P->module== "SELECTION"){
				if (P->p3.find(F) !=P->p3.end()){
					P->p3[F] 	= "1";
					P->N+=1;
				}else{
					F 			= "";
				}
			}else if (P->module== "BIDIR"){
				if (P->p4.find(F) !=P->p4.end()){
					P->p4[F] 	= "1";
					P->N+=1;
				}else{
					F 			= "";
				}
			}
		}
		else if (not F.empty()) {
			if (P->module=="MODEL"){
				if (P->p.find(F) !=P->p.end()){
					P->p[F]=string(*argv);
				}
			}else if(P->module=="FORMAT"){
				if (P->p2.find(F) !=P->p2.end()){
					P->p2[F]=string(*argv);
				}
			}else if(P->module=="SELECTION"){
				if (P->p3.find(F) !=P->p3.end()){
					
					P->p3[F]=string(*argv);
				}
			}else if(P->module=="BIDIR"){
				if (P->p4.find(F) !=P->p4.end()){
					
					P->p4[F]=string(*argv);
				}
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


