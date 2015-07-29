#include "model_selection.h"
#include "dirent.h"
#include <string>
#include <iostream>
#include "split.h"
#include <fstream>
#include <map>
#include "model.h"
#include "across_segments.h"
#include "template_matching.h"
using namespace std;
model::model(int k, double LL){
	K=k;
	ll=LL;
}
model::model(){}

S::S(string chr, int st, int sp, double n){
	chrom 	= chr;
	start 	= st;
	stop 	= sp;
	N 		= n;
}

void S::add_model(int k, double ll ){
	models[k].push_back(model(k, ll));
}

double BIC(double K, double ll, double N, double penality){
	return -2*ll + penality*(K+1)*6*log(N);
}

string S::print_best(){
	string text = "";
	for (int i = 0; i < best.bidirs.size(); i++){
		if (best.bidirs[i].w > 0){
			text+=("N: " + to_string(best.bidirs[i].mu) + "," + to_string(best.bidirs[i].si) + ","+
				to_string(best.bidirs[i].l) + "," + 
				to_string(best.bidirs[i].w) + "," + to_string(best.bidirs[i].pi) + "\n");
		}
	}
	for (int i = 0; i < best.elongations.size(); i++){
		if (best.elongations[i].w > 0.0000001){
			text+=("U: " + to_string(best.elongations[i].a) + "," + to_string(best.elongations[i].b) + ","+
				to_string(best.elongations[i].w) + "," + 
				to_string(best.elongations[i].pi)  + "\n");
		}
	}
	return text;
}

string S::print_header(){
	string text	= "#" + chrom + ":" + to_string(start) + "-" + to_string(stop)+"\n";
	return text;
}

void S::pick_best(double penality){
	typedef map<int , vector<model> >::iterator it_type;
	typedef map<int , model >::iterator it_type_2;
	double ll = nINF;
	for (it_type i = models.begin(); i != models.end(); i++){
		//pick best likelihood
		ll 		=nINF;
		for (int r = 0; r < i->second.size(); r++){
			if (i->second[r].ll > ll){
				ll 	= i->second[r].ll;
				bests[i->first] 	= i->second[r];
			}
		}
	}
	double score 	= INF;
	print_out 		= false;

	for (it_type_2 i = bests.begin(); i!= bests.end(); i++){
		if (BIC(i->first, i->second.ll, N, penality)  < score){
			score 		= BIC(i->first, i->second.ll, N, penality);
			best 		= i->second;
			print_out 	= true;
		}
	}

}

string S::getBestComponent(double x){
	double MAX 	= nINF;
	double prob;
	string out 	= "";
	for (int i = 0; i < best.bidirs.size(); i++){
		prob 	= max(best.bidirs[i].pdf(x, 1) , best.bidirs[i].pdf(x, -1) );
		if (prob > MAX ){
			MAX  	= prob;
			out 	= "bidir";
		}
	}
	for (int i = 0; i < best.elongations.size(); i++){
		prob 	= max(best.elongations[i].pdf(x, 1) , best.elongations[i].pdf(x, -1) );
		// if (prob and best.bidirs.size() > 0){
		// 	printf("%d,%f,%f \n", best.elongations[i].pi ,best.elongations[i].pdf(x, 1),best.elongations[i].pdf(x, -1));
		// }
		if (prob > MAX ){
			MAX  	= prob;
			if (best.elongations[i].pi > 0.5 and best.elongations[i].pos ==0 ){
				out 	= "forward";
			}else if (best.elongations[i].pi < 0.5 and best.elongations[i].pos ==0) {
				out 	= "reverse";
			}else{
				out 	= "noise";
			}
		}
	}
	return out;


}

void model::add_component(bool bidir , vector<string> line_array, int K){
	if (bidir){
		bidirs.push_back(EMG(stod(line_array[0]), stod(line_array[1]),
		stod(line_array[2]),stod(line_array[3]),stod(line_array[4]) ) );
	}else{
		if (K >0){
			elongations.push_back(UNI(stod(line_array[0]), stod(line_array[1]),
				stod(line_array[2]),stoi(line_array[3]), 0 ));
		}else{
			elongations.push_back(UNI(stod(line_array[0]), stod(line_array[1]),
				stod(line_array[2]),stoi(line_array[3]), 1));
		
		}
	}
}






map<string, vector<S> > read_in_file(string filename, map<string, vector<S> > model_fits){
	ifstream FH(filename);
	string line, chrom, stsp;
	string p 	 		= "#";
	string t 	 		= "~";
	string NN 	 	 	= "N: ";
	string UU 	 		= "U: ";

	const char *N_ID 	=NN.c_str();
	const char *U_ID 	=UU.c_str();
	
	
	const char *pound 	=p.c_str();
	const char *tilda 	=t.c_str();
	bool header 		= 1;
	double N, LL;
	vector<string> line_array ;
	int start, stop, K;
	while (getline(FH, line)){
		if (not header){
			if (line[0]==*pound){
				line_array 	= splitter(line.substr(1), ",");
				N 			= stod(line_array[1]);
				line_array 	= splitter(line_array[0], ":");
				chrom 		= line_array[0], stsp=line_array[1];
				line_array 	= splitter(stsp, "-");
				start=stoi(line_array[0]), stop=stoi(line_array[1]);
				model_fits[chrom].push_back(S(chrom, start, stop, N));
			}else if(line[0]==*tilda){
				line_array 	= splitter(line.substr(1), ",");
				K=stoi(line_array[0]), LL=stod(line_array[1]);
				model_fits[chrom].back().add_model(K,LL);
			}else if(line.substr(0,3) == NN or line.substr(0,3) == UU ){
				line_array 	= splitter(line.substr(3), ",");
				model_fits[chrom].back().models[K].back().add_component(line.substr(0,3) == NN, line_array, K);
			}
		}else{
			header=0;
		}
	}
	return model_fits;
	
	

}


bool checkIfFILE(string filename){
	ifstream FH(filename);
	string line;
	string pound 		= "#";
	const char *h 	=pound.c_str();
	if (FH) {
	  	getline(FH, line);
	  	getline(FH, line);
	  	if (line[0]==*h){
			FH.close();
	  		return true;
	  	}
	}
	FH.close();
	return false;
}

string check_out_files(string directory, string FILE, int i){
	string template_file 	= directory + to_string(i) + "_" + FILE ;
	ifstream FH(template_file);
	if (FH){
		FH.close();
		return check_out_files(directory, FILE, i+1);
	}
	FH.close();
	return template_file;
}


void run_model_selection(string in_directory, string out_directory, double penality){
	DIR *dir;
	struct dirent *ent;
	string FILE;
	map<string, vector<S> > model_fits ;
	if ((dir = opendir (in_directory.c_str())) != NULL) {
		/* print all the files and directories within directory */
		while ((ent = readdir (dir)) != NULL) {
			if (checkIfFILE(in_directory+string(ent->d_name))){
				FILE=string(ent->d_name);
				model_fits = read_in_file(in_directory+FILE, model_fits);
			}
		}
	}else {
		/* could not open directory */
		cout<<"couldn't open directory: "<<in_directory<<endl;
	}
	closedir (dir);
	
	string IGV_OUT 	= check_out_files(out_directory, "model_fits.bed", 1);
	string EMG_OUT 	= check_out_files(out_directory, "model_fits.txt", 1);
	
	ofstream FHW_IGV;
	ofstream FHW_EMG;
	
	FHW_IGV.open(IGV_OUT);
	FHW_EMG.open(EMG_OUT);
	FHW_IGV<<"track description=\"Mode Fits\"  cgColour1=white cgColour2=yellow cgColour3=red height=30\n";
	typedef map<string, vector<S> >::iterator it_type;
	//=============================
	//before model selection
	string 	current, prev_state, annotation, RGB,start, stop, ID, strand, INFO;
	int prev_start, dist;
	for (it_type i = model_fits.begin(); i!=model_fits.end(); i++){

		for (int j = 0 ; j < i->second.size();j++){
			i->second[j].pick_best(penality);
			FHW_EMG<<i->second[j].print_header();
			if (i->second[j].print_out){
				FHW_EMG<<i->second[j].print_best();
			}
			//now we want to output model fits to IGV
			
			if (i->second[j].print_out){
				for (int k = 0; k < i->second[j].best.elongations.size();k++){
					if (i->second[j].best.elongations[k].w > 0.0001){
						start 		= to_string(int(i->second[j].best.elongations[k].a*100 +
						 i->second[j].start));
						stop 		= to_string(int(i->second[j].best.elongations[k].b*100 +
						 i->second[j].start));
						if (i->second[j].best.elongations[k].pi == 1  ){
							RGB 	= "0,0,255";
							ID 		= "forward";
							strand 	= "+";
						}else{
							RGB 	= "255,0,0";
							ID 		= "reverse";
							strand 	= "-";
						}
						dist 		= int((i->second[j].best.elongations[k].b - i->second[j].best.elongations[k].a)*100);
						if (dist > 1000){
							annotation 	= (i->second[j].chrom + "\t" + start + "\t" + stop + "\t" + ID +
							"\t200\t" + strand + "\t"  + start + "\t" + stop + "\t" +  RGB +"\t" +"\t2\t100,100\t" + "0," + to_string(dist-100) +  "\n" );
							FHW_IGV<<annotation;	
				
						}
					}
				}

				for (int k = 0; k < i->second[j].best.bidirs.size();k++){
					if (i->second[j].best.bidirs[k].w > 0.0001){

						ID 			= "Bidir";
						RGB 		= "0,255,0";
						start 		= to_string(int((i->second[j].best.bidirs[k].mu - i->second[j].best.bidirs[k].si - (1.0 / i->second[j].best.bidirs[k].l))*100 + i->second[j].start   ));
						stop 		= to_string(int((i->second[j].best.bidirs[k].mu + i->second[j].best.bidirs[k].si + (1.0 / i->second[j].best.bidirs[k].l))*100 + i->second[j].start   ));
						INFO 		= to_string(int((i->second[j].best.bidirs[k].mu))) + "_" + to_string(int(i->second[j].best.bidirs[k].si*100)) + "_" + to_string(int(1.0/i->second[j].best.bidirs[k].l));
						annotation 	= (i->second[j].chrom + "\t" + start + "\t" + stop + "\t" + ID + "\t200\t.\t" + start + "\t" + stop+"\t" + RGB + "\t" + INFO + "\n");
						FHW_IGV<<annotation;	

						
					}
				}
			} 
		}
	}






	
}