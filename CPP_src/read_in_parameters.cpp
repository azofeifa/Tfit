#include <iostream>
#include <fstream>
#include <unistd.h>
#include <string>
#include <vector>
#include <algorithm>
#include "read_in_parameters.h"
#include <stdio.h>
#include <ctype.h>

#include <stdio.h>
#include <time.h>
#include "split.h"
using namespace std;

params::params(){
	p["-N"] 		= "EMG";
	p["-v"] 		= "1";
	p["-i"] 		= "";
	p["-j"] 		= "";
	p["-k"] 		= "";
	
	p["-o"] 		= "";
	p["-br"] 		= "300";
	p["-ns"] 		= "100";
	p["-minK"] 		= "1";
	p["-maxK"] 		= "3";
	p["-rounds"] 	= "10";
	p["-ct"] 		= "0.0001";
	p["-max_noise"] = "0.05";
	p["-chr"] 		= "all";
	p["-elon"] 		= "0";
	p["-mi"] 		= "300";
	p["-r_mu"] 		= "0";
	p["-ALPHA_0"] 	= "1";
	p["-BETA_0"] 	= "1";
	p["-ALPHA_1"] 	= "1";
	p["-BETA_1"] 	= "1";
	p["-ALPHA_2"] 	= "1";
	p["-ALPHA_3"] 	= "1";
	p["-template"] 	= "1";
	p["-bct"] 		= "0.95";
	p["-window_res"]= "10";
	p["-log_out"] 	= "";
	p["-ms_pen"] 	= "1";
	p["-topology"] 	= "1";

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
	p4["-N"] 				= "EMG";
	p4["-i"] 				= "";
	p4["-j"] 				= "";
	p4["-k"] 				= "";
	p4["-f"] 				= "";
	p4["-optimize"] 		= "0";
	p4["-o"] 				= "1";
	p4["-ns"] 				= "100";
	p4["-br"] 				= "50";
	p4["-density"] 			= "1000";
	p4["-window_res"] 		= "100";
	p4["-bct"] 				= "0.9";
	
	p4["-chr"] 				= "all";
	p4["-opt_res"] 			= "5";
	p4["-np"] 				= "4";
	p4["-MLE"] 				= "1";
	p4["-pad"] 				= "3000";

	p4["-rounds"] 			= "10";
	p4["-ct"] 				= "0.0001";
	p4["-max_noise"] 		= "0.05";
	p4["-mi"] 				= "300";
	p4["-r_mu"] 			= "0";
	p4["-ALPHA_0"] 			= "1";
	p4["-BETA_0"] 			= "1";
	p4["-ALPHA_1"] 			= "1";
	p4["-BETA_1"] 			= "1";
	p4["-ALPHA_2"] 			= "1";
	p4["-ALPHA_3"] 			= "1";
	p4["-elon"] 			= "0";
	p4["-show_seeds"] 		= "0";
	p4["-foot_res"] 		= "5";
	p4["-log_out"] 			= "";
	p4["-nf"] 				= "";


	p5["-v"] 			= "1";
	p5["-i"]  			= "";
	p5["-j"]  			= "";
	p5["-o"] 			= "";
	p5["-br"] 			= "300";
	p5["-ns"] 			= "100";
	p5["-minK"] 		= "1";
	p5["-maxK"] 		= "2";
	p5["-rounds"] 		= "10";
	p5["-r_mu"] 		= "1";
	p5["-mi"] 			= "1000";
	p5["-ct"] 			= "0.0001";
	p5["-max_noise"] 	= "0.05";
	p5["-np"] 			= "4";
	p5["-chr"] 			= "all";
	p5["-BETA_0"] 		= "1";
	p5["-ALPHA_0"] 		= "1";
	p5["-template"] 	= "1";
	p5["-opt_res"] 		= "10";
	p5["-show_seeds"] 	= "0";
	p5["-bct"] 			= "1";
	p5["-template"] 	= "0";
	p5["-pad"] 			= "0";
	
	p6["-v"] 			= "1";
	p6["-N"] 			= "";
	p6["-i"] 			= "";
	p6["-j"] 			= "";
	p6["-k"] 			= "";
	p6["-o"]			= "";
	p6["-chr"] 			= "all";
	p6["-br"] 			= "100";
	p6["-pad"] 			= "";
	p6["-ns"] 			= "100";
	p6["-r_mu"] 		= "1";
	p6["-ALPHA_0"] 		= "1";
	p6["-BETA_0"] 		= "1";
	p6["-ALPHA_1"] 		= "1";
	p6["-BETA_1"] 		= "1";
	p6["-ALPHA_2"] 		= "1";
	p6["-foot_res"] 	= "5";
	p6["-ALPHA_3"] 		= "1";
	p6["-rounds"] 		= "1";
	p6["-brounds"] 		= "1";
	p6["-mi"] 			= "1000";
	p6["-ct"] 			= "0.0001";
	p6["-max_noise"] 	= "0.05";
	p6["-log_out"] 		= ""	;

	
	N 				= 0;
	module 			= "";
	EXIT 			= 0;

}
const std::string currentDateTime() {
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
    // for more information about date/time format
    strftime(buf, sizeof(buf), "%m/%d/%Y %X", &tstruct);

    return buf;
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
void params::display(int nodes, int cores){
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
			cout<<"-r_mu     : "<<p["-r_mu"]<<endl;
			cout<<"-threads  : "<<cores<<endl;
			cout<<"-MPI_nodes: "<<nodes<<endl;
			cout<<"Questions/Bugs? joseph[dot]azofeifa[at]colorado[dot]edu"<<endl;
			cout<<"----------------------------------------------------------------"<<endl;
		}
		else{
			cout<<"----------------------------------------------------------------"<<endl;
			cout<<"              User Provided EMGU Parameters                     "<<endl;
			cout<<"                 (RUNNING MIXTURE MODEL)                        "<<endl;
			cout<<"            ...coupled to template matching...                  "<<endl;	
			cout<<"-i         : "<<p["-i"]<<endl;
			cout<<"-j         : "<<p["-j"]<<endl;
			cout<<"-k         : "<<p["-k"]<<endl;
			cout<<"-o         : "<<p["-o"]<<endl;
			cout<<"-log_out   : "<<p["-log_out"]<<endl;
			cout<<"-chr       : "<<p["-chr"]<<endl;
			cout<<"-br        : "<<p["-br"]<<endl;
			cout<<"-ns        : "<<p["-ns"]<<endl;
			cout<<"-rounds    : "<<p["-rounds"]<<endl;
			cout<<"-ct        : "<<p["-ct"]<<endl;
			cout<<"-mi        : "<<p["-mi"]<<endl;
			cout<<"-bct       : "<<p["-bct"]<<endl;
			cout<<"-window_res: "<<p["-window_res"]<<endl;
			cout<<"-minK      : "<<p["-minK"]<<endl;
			cout<<"-maxK      : "<<p["-maxK"]<<endl;			
			cout<<"-threads   : "<<cores<<endl;
			cout<<"-MPI_nodes : "<<nodes<<endl;
			cout<<"Questions/Bugs? joseph[dot]azofeifa[at]colorado[dot]edu"<<endl;
			cout<<"----------------------------------------------------------------"<<endl;
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
		cout<<"Questions/Bugs? joseph[dot]azofeifa[at]colorado[dot]edu"<<endl;
		cout<<"----------------------------------------------------------------"<<endl;
		
	}else if (module=="SELECTION"){
		cout<<"----------------------------------------------------------------"<<endl;
		cout<<"              User Provided EMGU Parameters                     "<<endl;
		cout<<"                    (MODEL SELECTION)                        "<<endl;
		cout<<"-i           : "<<p3["-i"]<<endl;
		cout<<"-o           : "<<p3["-i"]<<endl;
		cout<<"-penality    : "<<p3["-penality"]<<endl;
		cout<<"-to_igv      : "<<p3["-to_igv"]<<endl;
		cout<<"-to_EMG      : "<<p3["-to_EMG"]<<endl;
		cout<<"Questions/Bugs? joseph[dot]azofeifa[at]colorado[dot]edu"<<endl;
		cout<<"----------------------------------------------------------------"<<endl;
	}else if (module=="BIDIR")	{
		if (p4["-optimize"]=="0"){
			cout<<"----------------------------------------------------------------"<<endl;
			cout<<"              User Provided EMGU Parameters                     "<<endl;
			cout<<"                 (BIDIRECTIONAL DETECTOR)                       "<<endl;
			cout<<"-N           : "<<p4["-N"]<<endl;
			cout<<"-i           : "<<p4["-i"]<<endl;
			cout<<"-j           : "<<p4["-j"]<<endl;
			cout<<"-nf          : "<<p4["-nf"]<<endl;

			cout<<"-o           : "<<p4["-o"]<<endl;
			cout<<"-log_out     : "<<p4["-log_out"]<<endl;
			cout<<"-chr         : "<<p4["-chr"]<<endl;
			cout<<"-window_res  : "<<p4["-window_res"]<<endl;
			cout<<"-bct         : "<<p4["-bct"]<<endl;
			cout<<"-MLE         : "<<p4["-MLE"]<<endl;
			cout<<"-mi          : "<<p4["-mi"]<<endl;
			cout<<"-ct          : "<<p4["-ct"]<<endl;
			cout<<"-foot_res    : "<<p4["-foot_res"]<<endl;
			cout<<"-pad         : "<<p4["-pad"]<<endl;
			cout<<"-ns          : "<<p4["-ns"]<<endl;
			cout<<"-br          : "<<p4["-br"]<<endl;
			cout<<"-rounds      : "<<p4["-rounds"]<<endl;
			cout<<"-threads     : "<<cores<<endl;
			cout<<"-MPI_nodes   : "<<nodes<<endl;
			cout<<"Questions/Bugs? joseph[dot]azofeifa[at]colorado[dot]edu"<<endl;
			cout<<"----------------------------------------------------------------"<<endl;
		
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
			cout<<"-ns          : "<<p4["-ns"]<<endl;
			cout<<"-br          : "<<p4["-br"]<<endl;
			cout<<"-threads     : "<<cores<<endl;
			cout<<"-MPI_nodes   : "<<nodes<<endl;
			cout<<"Questions/Bugs? joseph[dot]azofeifa[at]colorado[dot]edu"<<endl;
			cout<<"----------------------------------------------------------------"<<endl;
			
		}
	}else if (module == "SINGLE") {
		cout<<"----------------------------------------------------------------"<<endl;
		cout<<"              User Provided EMGU Parameters                     "<<endl;
		cout<<"                 (Fitting single model)                       "<<endl;
		cout<<"Date Time    : "<<currentDateTime()<<endl;
		cout<<"-i           : "<<p5["-i"]<<endl;
		cout<<"-j           : "<<p5["-j"]<<endl;
		cout<<"-o           : "<<p5["-o"]<<endl;
		cout<<"-ns          : "<<p5["-ns"]<<endl;
		cout<<"-br          : "<<p5["-br"]<<endl;
		cout<<"-bct         : "<<p5["-bct"]<<endl;
		cout<<"-opt_res     : "<<p5["-opt_res"]<<endl;
		cout<<"-chr         : "<<p5["-chr"]<<endl;
		cout<<"-template    : "<<p5["-template"]<<endl;
		cout<<"-threads     : "<<cores<<endl;
		cout<<"-MPI_nodes   : "<<nodes<<endl;
		cout<<"Questions/Bugs? joseph[dot]azofeifa[at]colorado[dot]edu"<<endl;
		cout<<"----------------------------------------------------------------"<<endl;
	}else if (module == "BOOTSTRAP"){

		cout<<"----------------------------------------------------------------"<<endl;
		cout<<"              User Provided EMGU Parameters                     "<<endl;
		cout<<"               (Bootstraping MLE Estimates)                     "<<endl;
		cout<<"Date Time    : "<<currentDateTime()<<endl;
		cout<<"-i           : "<<p6["-i"]<<endl;
		cout<<"-j           : "<<p6["-j"]<<endl;
		cout<<"-k           : "<<p6["-k"]<<endl;
		cout<<"-o           : "<<p6["-o"]<<endl;
		cout<<"-log_out     : "<<p6["-log_out"]<<endl;
		cout<<"-ns          : "<<p6["-ns"]<<endl;
		cout<<"-br          : "<<p6["-br"]<<endl;
		cout<<"-chr         : "<<p6["-chr"]<<endl;
		cout<<"-rounds      : "<<p6["-rounds"]<<endl;
		cout<<"-brounds     : "<<p6["-brounds"]<<endl;
		
		cout<<"-pad         : "<<p6["-pad"]	<<endl;
		cout<<"-r_mu        : "<<p6["-r_mu"]	<<endl;
		cout<<"-ALPHA_0     : "<<p6["-ALPHA_0"]	<<endl;
		cout<<"-BETA_0      : "<<p6["-BETA_0"]	<<endl;
		cout<<"-ALPHA_1     : "<<p6["-ALPHA_1"]	<<endl;
		cout<<"-BETA_1      : "<<p6["-BETA_1"]	<<endl;
		cout<<"-ALPHA_2     : "<<p6["-ALPHA_2"]	<<endl;
		cout<<"-ALPHA_3     : "<<p6["-ALPHA_3"]	<<endl;
		cout<<"-mi          : "<<p6["-mi"]	<<endl;
		cout<<"-ct          : "<<p6["-ct"]	<<endl;
		cout<<"-max_noise   : "<<p6["-max_noise"]	<<endl;
		cout<<"-threads     : "<<cores<<endl;
		cout<<"-MPI_nodes   : "<<nodes<<endl;
		cout<<"Questions/Bugs? joseph[dot]azofeifa[at]colorado[dot]edu"<<endl;
		cout<<"----------------------------------------------------------------"<<endl;

	}
}


string params::get_header(int ID){
	string header = "";


	if (ID == 5){
		header+="#Date Time    : "+currentDateTime()+"\n";
		header+="#-i           : "+p5["-i"]+"\n";
		header+="#-j           : "+p5["-j"]+"\n";
		header+="#-o           : "+p5["-o"]+"\n";
		header+="#-ns          : "+p5["-ns"]+"\n";
		header+="#-br          : "+p5["-br"]+"\n";
		header+="#-chr         : "+p5["-chr"]+"\n";
		header+="#-mi          : "+p5["-mi"]+"\n";
		header+="#-ct          : "+p5["-ct"]+"\n";
		header+="#-rounds      : "+p5["-rounds"]+"\n";
		header+="#-pad         : "+p5["-pad"]+"\n";	
	}else if (ID == 4){
		header+="#----------------------------------------------------\n";
		header+="#Date Time    : "+currentDateTime()+"\n";
		header+="#-N           : "+p4["-N"]+"\n";
		header+="#-i           : "+p4["-i"]+"\n";
		header+="#-j           : "+p4["-j"]+"\n";
		header+="#-nf          : "+p4["-nf"]+"\n";
		header+="#-o           : "+p4["-o"]+"\n";
		header+="#-ns          : "+p4["-ns"]+"\n";
		header+="#-br          : "+p4["-br"]+"\n";
		header+="#-bct         : "+p4["-bct"]+"\n";
		header+="#-window_res  : "+p4["-window_res"]+"\n";
		header+="#-mi          : "+p4["-mi"]+"\n";
		header+="#-ct          : "+p4["-ct"]+"\n";
		header+="#-rounds      : "+p4["-rounds"]+"\n";
		header+="#-foot_res    : "+p4["-foot_res"]+"\n";
		header+="#-pad         : "+p4["-pad"]+"\n";	
	
		header+="#-ALPHA_0     : "+p4["-ALPHA_0"]+"\n";	
		header+="#-BETA_0      : "+p4["-BETA_0"]+"\n";	
		header+="#-BETA_1      : "+p4["-BETA_1"]+"\n";	
		header+="#-ALPHA_2     : "+p4["-ALPHA_2"]+"\n";	
		header+="#-ALPHA_3     : "+p4["-ALPHA_3"]+"\n";	

		header+="#----------------------------------------------------\n";
	}else if(ID==0){
		header+="#----------------------------------------------------\n";
		header+="#Date Time    : "+currentDateTime()+"\n";
		header+="#-N           : "+p["-N"]+"\n";
		header+="#-i           : "+p["-i"]+"\n";
		header+="#-j           : "+p["-j"]+"\n";
		header+="#-k           : "+p["-k"]+"\n";
		header+="#-o           : "+p["-o"]+"\n";
		header+="#-ns          : "+p["-ns"]+"\n";
		header+="#-br          : "+p["-br"]+"\n";
		header+="#-bct         : "+p["-bct"]+"\n";
		header+="#-window_res  : "+p["-window_res"]+"\n";
		header+="#-minK        : "+p["-minK"]+"\n";
		header+="#-maxK        : "+p["-maxK"]+"\n";
		header+="#-ms_pen      : "+p["-ms_pen"]+"\n";
		header+="#-mi          : "+p["-mi"]+"\n";
		header+="#-ct          : "+p["-ct"]+"\n";
		header+="#-rounds      : "+p["-rounds"]+"\n";
		header+="#-ALPHA_0     : "+p["-ALPHA_0"]+"\n";	
		header+="#-BETA_0      : "+p["-BETA_0"]+"\n";	
		header+="#-BETA_1      : "+p["-BETA_1"]+"\n";	
		header+="#-ALPHA_2     : "+p["-ALPHA_2"]+"\n";	
		header+="#-ALPHA_3     : "+p["-ALPHA_3"]+"\n";	
		header+="#----------------------------------------------------\n";
	
	}else if (ID==6){
		header+="#----------------------------------------------------\n";
		header+="#Date Time    : "+currentDateTime()+"\n";
		header+="#-N           : "+p6["-N"]+"\n";
		header+="#-i           : "+p6["-i"]+"\n";
		header+="#-j           : "+p6["-j"]+"\n";
		header+="#-k           : "+p6["-k"]+"\n";
		header+="#-o           : "+p6["-o"]+"\n";
		header+="#-ns          : "+p6["-ns"]+"\n";
		header+="#-br          : "+p6["-br"]+"\n";
		header+="#-mi          : "+p6["-mi"]+"\n";
		header+="#-ct          : "+p6["-ct"]+"\n";
		header+="#-foot_res    : "+p6["-foot_res"]+"\n";
		header+="#-rounds      : "+p6["-rounds"]+"\n";
		header+="#-brounds     : "+p6["-brounds"]+"\n";
		header+="#-ALPHA_0     : "+p6["-ALPHA_0"]+"\n";	
		header+="#-BETA_0      : "+p6["-BETA_0"]+"\n";	
		header+="#-BETA_1      : "+p6["-BETA_1"]+"\n";	
		header+="#-ALPHA_2     : "+p6["-ALPHA_2"]+"\n";	
		header+="#-ALPHA_3     : "+p6["-ALPHA_3"]+"\n";	
		header+="#----------------------------------------------------\n";

	}

	return header;
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
				}else if(P->module =="SINGLE"){
					if (P->p5.find(ID) !=P->p5.end()){
							P->p5[ID] 	= val;
					}else if (not ID.empty() ){
						not_valid.push_back(ID);
					}
				}else if(P->module=="BOOTSTRAP"){
					if (P->p6.find(ID) !=P->p6.end()){
							P->p6[ID] 	= val;
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
			}else if (P->module == "SINGLE"){
				if (P->p5.find(F) !=P->p5.end()){
					P->p5[F] 	= "1";
					P->N+=1;
				}else{
					F 			= "";
				}
			}else if (P->module == "BOOTSTRAP"){
				if (P->p6.find(F) !=P->p6.end()){
					P->p6[F] 	= "1";
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
			}else if(P->module=="SINGLE"){
				if (P->p5.find(F) !=P->p5.end()){
					P->p5[F]=string(*argv);
				}
			}else if(P->module=="BOOTSTRAP"){
				if (P->p6.find(F) !=P->p6.end()){
					P->p6[F]=string(*argv);
				}
			}
		}
		argv++;
		GO_FORIT=true;
	}

}

void fill_in_bidir_boostrap(params * P ){
	string FILE 		= P->p6["-i"];
	ifstream FH(FILE);
	vector<string> lineArray;
	string line, ID;
	map<string, int> NOT;
	NOT["-i"] 	= 1;
	NOT["-j"] 	= 1;
	NOT["-N"] 	= 1;
	NOT["-o"] 	= 1;
	NOT["-chr"] = 1;
	NOT["-log_out"] = 1;
	
	
	bool COLLECT;
	while(getline(FH, line)){
		ID 	= "";
		COLLECT 	= false;
		if (line.substr(0,1)=="#"){
			for (int i = 1; i < line.size(); i++){
				if (not isspace(line[i]) ){
					ID+=line[i];
				}else{
					break;
				}
			}
			if (P->p6.find(ID) != P->p6.end() and NOT.find(ID)==NOT.end() ){
				COLLECT=true;
			}
			if (COLLECT){
				lineArray 	= splitter(line, ":")	;
				if (lineArray.size()==2){
					P->p6[ID] 	= lineArray[1].substr(1,lineArray[1].size() );
				}
			}


		}else{

		}

	}
}


params * readInParameters( char* argv[]){	
	string userModParameter = "";
	params 	* P = new params;
	fillInOptions(argv, P);
	return P;
}


