import merge_data_types as mdt
import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.append("/Users/joeyazo/Desktop/Lab/EMG/python_src/")
import model
import simulate
def merge_data_out(FILE):
	G=list()
	with open(FILE) as FH:
		for line in FH:
			if "#" == line[0]:
				chrom,info 			= line[1:].strip("\n").split(":")
				start_stop, N,aN 	= info.split(",")
				start,stop 			= start_stop.split("-")
				G.append(mdt.segment(chrom,int(start),int(stop),float(N), annotation_N=int(aN)))
			elif "~" == line[0]:
				G[-1].insert_model_info(line)
			elif "N:"==line[:2] or "U:"==line[:2]:
				G[-1].insert_component(line)
			else:
				line_array 				= line.strip("\n").split(",")
				data_type,peak, data 	= line_array[0], line_array[1],",".join(line_array[2:])
				if data_type != "dbSNP":
					data 					= [(float(d.split(",")[0]),float(d.split(",")[1])) for d in data.split(":") ]
				else:
					data 					= [(float(d.split(",")[0]), d.split(",")) for d in data.split(":")  ]
				setattr(G[-1], data_type, data)
				setattr(G[-1], data_type+"_peak", bool(peak))
	return G

if __name__ == "__main__":
	FILE 	= "/Users/joeyazo/Desktop/Lab/gro_seq_files/HCT116/merged_data_file_100.txt"
	G 		= merge_data_out(FILE)
	