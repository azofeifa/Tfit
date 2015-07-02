#include "write.h"
#include "load.h"
#include <string>
#include <iostream>
#include <fstream>
#include <map>
using namespace std;

void write_out_file(string out_file_name, annotations * A, bool write_data){
	ofstream FHW;
	FHW.open(out_file_name);
	map<string, annotation_cluster *> collections 	= A->collections;
	typedef map<string, annotation_cluster *>::iterator it_type;
	annotation_cluster * cluster;
	annotation * current;
	
	for(it_type iterator 	= collections.begin(); iterator != collections.end(); iterator++) {
		cluster 			= iterator->second;
		for (cluster = iterator->second; cluster!=NULL; cluster=cluster->next){
			for (current = cluster->root; current!=NULL; current=current->next){
				if (current->rv_dgu_best != NULL){
					FHW<<"#"<<current->name<<","<<current->chrom<<":"<<to_string(current->start)<<"-"<<to_string(current->stop)<<":"<<current->strand<<","<<to_string(current->density)<<","<<to_string(current->total)<<endl;
					FHW<<to_string(current->rv_dgu_best->i)<<","<<to_string(current->rv_dgu_best->u)<<","<<to_string(current->rv_dgu_best->s)<<","<<to_string(current->rv_dgu_best->l)<<","<<to_string(current->rv_dgu_best->w)<<endl;
					if (write_data){
						for (int i =0; i < current->N;i++){
							FHW<<to_string(current->D[i][0])<<","<<to_string(current->D[i][1])<<","<<to_string(current->D[i][2])<<","<<endl;
						}
					}
				}
			}
		}
	}


}