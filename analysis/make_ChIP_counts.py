import sys, BIC, load
import math
import numpy as np
import matplotlib.pyplot as plt
def check_bidir_component(rv, si_thresh=4, l_thresh=2., w_thresh=0.01, pi_thresh=0.1):
	if rv.type == "N":
		if rv.si < si_thresh and rv.l < l_thresh and rv.w > w_thresh and pi_thresh <= rv.pi<= (1-pi_thresh):
			return True
	return False
def bin_ChIP_signal(N, I,res=100, limit=4):
	
	std 		= N.si + math.sqrt(pow(1.0 / N.l, 2))
	bins 		= np.linspace(-4,4,100)*std
	bins 		+=N.mu
	XS 			= list()
	for data_type in I.data_types:
		X 		= np.zeros((len(bins), 2))
		X[:,0] 	= bins
		j 		= 0
		NN 		= len(bins)
		getattr(I, data_type).sort()
		for x,y in getattr(I, data_type):
			while j < NN and X[j,0] <= x:
				j+=1
			if data_type!="dbSNP" and 1< j < NN and x < X[j,0]:
				X[j-1, 1]+=y
			elif 1< j < NN and  x < X[j,0]:
				X[j-1, 1]+=1
		XS.append((X,data_type,getattr(I, data_type + "_peak")))
	return XS

def run(G, penality, diff_threshold, out_file_name, si_thresh, l_thresh, w_thresh, pi_thresh):
	FHW 	= open(out_file_name+"_" + str(penality) + "_" + str(diff_threshold) +"_" + str(si_thresh) + "_" + str(l_thresh) + "_" + str(w_thresh)+ "_" + str(pi_thresh)   ,"w"  )
	for I in G:
		model 	= BIC.get_best_model(I, penality , diff_threshold)
		bidirs 	= [rv for rv in model.rvs if check_bidir_component(rv,si_thresh=si_thresh, l_thresh=l_thresh, w_thresh=w_thresh, pi_thresh=pi_thresh)]
		if bidirs:
			FHW.write("#" + I.chrom + ":" + str(I.start) + "-" + str(I.stop)+ "\n")

			for N in bidirs:
				XS 	= bin_ChIP_signal(N, I)
				FHW.write(N.__str__() + "\n")
				for X, data_type, peak in XS:
					FHW.write(data_type+"," + str(peak) + "," + ",".join([str(x) +"-"+str(y) for x,y in X] ) + "\n")







if __name__ == "__main__":
	if len(sys.argv)==1:
		merge_data_out 		= "/Users/joeyazo/Desktop/Lab/gro_seq_files/HCT116/merged_data_file_20.txt"
		out_file_name 		= "/Users/joeyazo/Desktop/Lab/gro_seq_files/HCT116/counted_data"
		G 					= load.merge_data_out(merge_data_out)
		penality 			= 100
		diff_threshold 		= 5
		si_thresh, l_thresh, w_thresh, pi_thresh 	= 4,2,0.01, 0.1

		run(G, penality, diff_threshold, out_file_name,si_thresh, l_thresh, w_thresh, pi_thresh)
