import time, load, merge_data_types as mdt
import sys, BIC
import matplotlib.pyplot as plt
import numpy as np
def output(I,  FHW, penality,diff_threshold ):
	model 	= BIC.get_best_model(I, penality , diff_threshold)
	FHW.write("#" + I.chrom + ":" + str(I.start) + "-" + str(I.stop) +  "," + str(I.annotation_N) + "\n")
	for rv in model.rvs:
		FHW.write(rv.__str__()+"\n")

def run(merged_file, out_file_name, penality,diff_threshold):
	FHW 	= open(out_file_name+"_" + str(penality) + "_" + str(diff_threshold) ,"w"  )
	I 		= None
	with open(merged_file) as FH:
		for line in FH:
			if "#" == line[0]:
				if I is not None:
					output(I, FHW, penality,diff_threshold)
				chrom,info 			= line[1:].strip("\n").split(":")
				start_stop, N,aN 	= info.split(",")
				start,stop 			= start_stop.split("-")
				I 					= mdt.segment(chrom,int(start),int(stop),float(N), annotation_N=int(aN))
			elif "~" == line[0]:
				I.insert_model_info(line)
			elif "N:"==line[:2] or "U:"==line[:2]:
				I.insert_component(line)
	FHW.close()
def read_in_display(FILE, w=0.1, si=True, NN=1):
	with open(FILE) as FH:
		collect=False
		G 		= {"mu":list(), "si":list(), "l":list(), "w":list(), "pi":list()}
		TYPES 	= ("mu","si", "l", "w", "pi" )
		for line in FH:
			if "#" == line[0]:
				N 	= int(line.strip("\n").split(",")[-1])
				collect=False
				if not si or N==NN:
					collect=True
			elif collect and line[:2] == "N:":
				ps 					= [float(p) for p in line[3:].strip("\n").split(",")]
				if ps[2] != 2.:
					for p, ty in zip(ps, TYPES):
						G[ty].append(p)

	si_mean 	= np.mean([s*100 for s in G["si"] if s< 10])
	counts,edges 	= np.histogram([s*100 for s in G["si"] if s< 10], bins=200)
	edges 		= (edges[:-1] + edges[1:]) / 2.
	si_mode 	= max([(x,y) for x,y in zip(counts, edges)])[1]
	l_mean 		= np.mean([100. / (l) for l in G["l"] if 100. / (l) < 1000])
	counts,edges 	= np.histogram([100. / (l) for l in G["l"] if 100. / (l) < 1000], bins=200)
	edges 		= (edges[:-1] + edges[1:]) / 2.
	l_mode 		= max([(x,y) for x,y in zip(counts, edges)])[1]
	pi_mean 	= np.mean(G["pi"])
	
	return si_mean, si_mode, l_mean, l_mode, pi_mean, len(G["pi"])


if __name__ == "__main__":
	RUN 				= False
	if RUN:
		if len(sys.argv)==1:
			merged_file 	= "/Users/joeyazo/Desktop/Lab/gro_seq_files/HCT116/merged_data_file_100.txt"
			out 			= "/Users/joeyazo/Desktop/BIC_BEST"
			penality 		= 100
			diff_threshold 	= 10
		else:
			
			merged_file 	= sys.argv[1]
			out 			= sys.argv[2]
			penality 		= float(sys.argv[3])
			diff_threshold 	= float(sys.argv[4])
		
		run(merged_file, out, penality,diff_threshold)
	else:
		FILE 			= "/Users/joeyazo/Desktop/Lab/EMG_files/final_models_100.0_3.0"
		print read_in_display(FILE , w=0.1, si=True, NN=1)
		print read_in_display(FILE , w=0.1, si=True, NN=0)
		

