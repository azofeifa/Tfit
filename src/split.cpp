#include "split.h"
#include <iostream>
#include <sstream>
using namespace std;






vector<string> string_split(string s, const char delimiter)
{
  size_t start=0;
  size_t end=s.find_first_of(delimiter);
  
  vector<string> output;
  
  while (end <= string::npos)
    {
      output.emplace_back(s.substr(start, end-start));
      
      if (end == string::npos)
	break;
      
      start=end+1;
      end = s.find_first_of(delimiter, start);
    }
  
  return output;
}



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
vector<string> splitter2(string line, string delim){
	vector<string> tokens;
	istringstream iss(line);
	string token;
	while(std::getline(iss, token, '\t' )){   // but we can specify a different one
		tokens.push_back(token);
	}
	return tokens;
}

vector<string> split_by_bar(string line, string delim){
	vector<string> tokens;
	istringstream iss(line);
	string token;
	while(std::getline(iss, token, '|' )){   // but we can specify a different one
		tokens.push_back(token);
	}
	return tokens;
}
vector<string> split_by_colon(string line, string delim){
	vector<string> tokens;
	istringstream iss(line);
	string token;
	while(std::getline(iss, token, ':' )){   // but we can specify a different one
		tokens.push_back(token);
	}
	return tokens;
}
vector<string> split_by_tab(string line, string delim){
	vector<string> tokens;
	istringstream iss(line);
	string token;
	while(std::getline(iss, token, '\t' )){   // but we can specify a different one
		tokens.push_back(token);
	}
	return tokens;
}
vector<string> split_by_comma(string line, string delim){
	vector<string> tokens;
	istringstream iss(line);
	string token;
	while(std::getline(iss, token, ',' )){   // but we can specify a different one
		tokens.push_back(token);
	}
	return tokens;
}
vector<string> split_by_dash(string line, string delim){
	vector<string> tokens;
	istringstream iss(line);
	string token;
	while(std::getline(iss, token, '-' )){   // but we can specify a different one
		tokens.push_back(token);
	}
	return tokens;
}


string strip(string ELE, string D){
	const char *d 	= D.c_str();
	string result 	= "";
	for (int i = 0; i < ELE.size(); i++){
		if (ELE[i]==*d){
			break;
		}else{
			result+=ELE[i];
		}
	}
	return result;
}

string join(vector<string> toBeJoined, string delim){
	typedef vector<string>::iterator vs_it;
	string result="";
	for (vs_it i =toBeJoined.begin(); i != toBeJoined.end(); i++){
		result=result + delim + *i;
	}
	return result;
}	
