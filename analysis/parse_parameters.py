
def write_out(G, out, penality,diff_threshold ):
	FHW 	= open(out+"_" + str(penality) + "_" +str(diff_threshold) +".txt" , "w")
	for I in G:
		model 	= BIC.get_best_model(I, penality , diff_threshold)
		FHW.write("#" + I.chrom + ":" + str(I.start) + "-" + str(I.stop) + "\n")
		for rv in model.rvs:
			FHW.write(rv.__str__()+"\n")
	FHW.close()




if __name__ == "__main__":
	import sys
	if len(sys.argv)==1:
		merged_file 	= "/Users/joeyazo/Desktop/Lab/gro_seq_files/HCT116/merged_data_file_100.txt"
		out 			= "/Users/joeyazo/Desktop/BIC_BEST"
		penality 		= 100
		diff_threshold 	= 10
		EMG_path 		= "/Users/joeyazo/Desktop/Lab/EMG/python_src/"
	else:
		
		merged_file 	= sys.argv[1]
		out 			= sys.argv[2]
		penality 		= sys.argv[3]
		diff_threshold 	= sys.argv[4]
		EMG_path 		= sys.argv[5]
	
	sys.path.append(EMG_path)
	import time, load, merge_data_types as mdt
	import sys, BIC

	G 				= load.merge_data_out(merged_file)
	write_out(G, out, penality,diff_threshold)

