#ifndef split_H
#define split_H

#include <vector>
#include <string>
using namespace std;

vector<string> splitter(string ELE, string D);
string strip(string ELE, string D);
string join(vector<string>, string);
vector<string> splitter2(string ELE, string D);
vector<string> split_by_bar(string ELE, string D);
vector<string> split_by_colon(string ELE, string D);
vector<string> split_by_tab(string ELE, string D);
vector<string> split_by_comma(string ELE, string D);
vector<string> split_by_dash(string ELE, string D);


#endif