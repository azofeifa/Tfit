import sys
sys.path.append("/Users/azofeifa/Lab/")
sys.path.append("/Users/joazofeifa/Lab/")
from interval_searcher import intervals, node
import numpy as np
import matplotlib.pyplot as plt

def get_TSS_tree(FILE):
	G 	= {}
	with open(FILE) as FH:
		for line in FH:
			chrom, start, stop 	= line.strip("\n").split("\t")[:3]
			start, stop 		= int(start), int(stop)
			if chrom not in G:
				G[chrom] 	= list()
			G[chrom].append((start, stop))
	for chrom in G:
		G[chrom].sort()
		G[chrom]=node.tree(G[chrom])
	return G

def load(FILE):
	G ={}
	with open(FILE) as FH:
		for line in FH:
			if "#"!=line[0]:
				chrom,start, stop,ps 	= line.strip("\n").split("\t")
				if chrom not in G:
					G[chrom]=list()
				w 	= float(ps.split("_")[3])
				if w > 0.3	:
					G[chrom].append((int(start), int(stop), 1.0-w  ))
	for chrom in G:
		G[chrom].sort()
	return G
def move_elongation(B, TSS_T):
	weights 	= np.linspace(0.3,0.6,30)
	D 			= list()
	for weight in weights:
		N 	= 0.0
		F 	= 0.0
		for chrom in B:
			T 	= TSS_T[chrom]
			for start, stop, w in B[chrom]:
				if w > weight   :
					N+=1
					FINDS 	= T.searchInterval((start, stop))
					if FINDS:
						F+=1
		D.append((N,F))
	return D,weights
def get_box(B, TSS_T):
	genes, eRNAS 	= list(), list() 
	for chrom in B:
		T 	= TSS_T[chrom]
		for start, stop, w in B[chrom]:
			FINDS 	= T.searchInterval((start, stop))
			if FINDS:
				genes.append(w)
			else:
				eRNAS.append(w)
	plt.hist(genes,bins=30,alpha=0.5,label="genes",normed=1)
	plt.hist(eRNAS,bins=30,alpha=0.5,label="eRNAs",normed=1)
	plt.legend()
	plt.show()



def display_accuracy(D, weights):
	F 	= plt.figure(figsize=(10,7))
	ax 	= F.add_subplot(1,1,1)
	ax.scatter(weights, [F/N for N,F in D])
	plt.show()




if __name__ == "__main__":
	FILE 	= "/Users/joazofeifa/Lab/gro_seq_files/Allen2014/EMG_out_files/Allen2014_DMSO2_3-4_bidirectional_hits_intervals.bed"
	TSS_T 	= get_TSS_tree("/Users/joazofeifa/Lab/genome_files/TSS.bed")
	B 		= load(FILE)
	#D,weights 		= move_elongation(B, TSS_T)
	get_box(B,TSS_T)
	#display_accuracy(D,weights)
