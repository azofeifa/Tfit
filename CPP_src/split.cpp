#include "split.h"
#include <iostream>
using namespace std;

vector<string> splitter(string ELE, string D){
	int j = 0;
	vector<string> results;
	const char *d =D.c_str();
	while (not ELE.empty() and j != ELE.size()){
		if (ELE[j] == *d){
			results.push_back(ELE.substr(0,j));
			ELE=ELE.substr(j+1,ELE.size());
			j=0;
		}
		j++;
	}
	results.push_back(ELE.substr(0,j));

	return results;
}
string join(vector<string> toBeJoined, string delim){
	typedef vector<string>::iterator vs_it;
	string result="";
	for (vs_it i =toBeJoined.begin(); i != toBeJoined.end(); i++){
		result=result + delim + *i;
	}
	return result;
}	
