#ifndef error_stdo_logging_H
#define error_stdo_logging_H
#include <fstream>
#include <string>
#include <iostream>
using namespace std;
class Log_File{
public:
	int job_ID;
	string job_name;
	int rank; 
	string log_out_dir;
	ofstream FHW;
	Log_File();
	Log_File(int, int, string, string);
	void write(string, int);
};


#endif
