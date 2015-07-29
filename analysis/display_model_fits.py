import merge_data_types as mdt
import BIC
import matplotlib.pyplot as plt
import numpy as np
#=============================================
#GLOBAL VARIABLES
si_thresh 	= 4
l_thresh 	= 2
w_thresh 	= 0.01
pi_thresh 	= 0.1
#=============================================
penality 		= 75
diff_threshold  = 5

def output(I, forward, reverse):
	D 	= (list(), list())
	model 	= BIC.get_best_model(I, penality , diff_threshold)
	print model.k
	for rv in model.rvs:
		print rv
	for i,FH in enumerate((forward, reverse)):
		for line in FH:
			chrom,start, stop, cov 	= line.strip("\n").split("\t")
			pos 					= (float(stop) + float(start ) ) /2.
			if pos > I.stop:
				break
			elif I.start<=pos <=I.stop:

				D[i].append((pos, float(cov)))
	minX,maxX 		= min([x for d in D for x,y in d]), max([x for d in D for x,y in d])

	xs 				= np.linspace(0, (maxX-minX)/100., 1000 )
	bins 			= 500
	counts,edges 	= np.histogram([(x-minX)/100. for x,y in D[0]], weights=[y for x,y in D[0]], bins=bins, normed=1)
	edges 			= (edges[:-1] + edges[1:])/2.
	F 				= plt.figure(figsize=(15,10))
	plt.bar(edges, counts,width=(edges[-1]-edges[0])/bins,alpha=0.5)
	plt.plot(xs, map(lambda x: model.pdf(x, 1) , xs), linewidth=2.)

	counts,edges 	= np.histogram([(x-minX)/100. for x,y in D[1]], weights=[y for x,y in D[1]], bins=bins, normed=1)
	edges 			= (edges[:-1] + edges[1:])/2.
	plt.bar(edges, -counts, width=(edges[-1]-edges[0])/bins, color="red",alpha=0.5)
	plt.plot(xs, map(lambda x: -model.pdf(x, -1) ,xs ), linewidth=2.)
	
	plt.show()
	


def parse_file(FILE , forward, reverse):
	FH1, FH2 	= open(forward), open(reverse)
	with open(FILE) as FH:
		header,I	= True,None
		for line in FH:
			if header:
				if "#" == line[0]:
					if I is not None:
						output(I, FH1, FH2 )
					chrom,info 			= line[1:].strip("\n").split(":")
					start_stop, N  		= info.split(",")
					start,stop 			= start_stop.split("-")
					I 					= mdt.segment(chrom,int(start),int(stop),float(N) )
				elif "~" == line[0]:
					I.insert_model_info(line)
				elif "N:"==line[:2] or "U:"==line[:2]:
					I.insert_component(line)
			else:
				header=False


if __name__ == "__main__":
	OLD 	= False

	if OLD:
		EMG_IN 	= "/Users/joeyazo/Desktop/Lab/gro_seq_files/HCT116/EMG_out_files/test_file_2.tsv"
		EMG_OUT = "/Users/joeyazo/Desktop/Lab/gro_seq_files/HCT116/EMG_out_files/model_fits_out_all_11"

		run(EMG_IN, EMG_OUT)
	else:
		forward 	= "/Users/joeyazo/Desktop/Lab/gro_seq_files/HCT116/bed_graph_files/DMSO2_3.pos.BedGraph"
		reverse 	= "/Users/joeyazo/Desktop/Lab/gro_seq_files/HCT116/bed_graph_files/DMSO2_3.neg.BedGraph"

		fits 		= "/Users/joeyazo/Desktop/Lab/gro_seq_files/HCT116/EMG_out_files/DMSO_ND_intervals_model_fits_2/model_fits_out_chr1_1"
		parse_file(fits, forward, reverse)
	pass





