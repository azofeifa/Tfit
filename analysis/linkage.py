import time, sys
sys.path.append("/Users/joazofeifa/Lab/interval_searcher/")
sys.path.append("/Users/joazofeifa/Lab/")
from interval_searcher import intervals, node
import numpy as np
import matplotlib.pyplot as plt
import math
def load(FILE, bidir_bootstrap=False, test=False,FILTER=None, isolate=None):
	G,t 	= {},0
	with open(FILE) as FH:
		for line in FH:
			if test and t > 100:
				break
			t+=1
			if line[0]!= "#":
				chrom,start, stop 	= line.split("\t")[:3]
				start, stop 		= int(start), int(stop)
				if FILTER is None and isolate is None or ( FILTER is not None and chrom in FILTER and not FILTER[chrom].searchInterval((int(start), int(stop)))) or (isolate is not None and chrom in isolate and isolate[chrom].searchInterval((int(start), int(stop)))) :
					if chrom not in G:
						G[chrom] 		= list()
					if not bidir_bootstrap:
						G[chrom].append((int(start), int(stop)))
					else:
						ci 				= line.split("\t")[4].split("_")[0]
						G[chrom].append((int(start) , int(stop) , float(ci)))					

	for chrom in G:
		G[chrom].sort()
	return G
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
def compute_distances(LST):
	j 	= 0
	D 	= list()
	while (j + 1) < len(LST):

		D.append(math.log(LST[j+1][1] - LST[j][0] ,10))
		j+=1
	return D
def compute_distances_random(LST):
	LST.sort()
	j 	= 0
	D 	= list()
	while (j + 1) < len(LST):

		D.append(math.log(LST[j+1] - LST[j]+1,10))
		j+=1
	return D
def linkage(LST,threshold=1000, single=False):
	LST.sort()
	j,N 	= 0,len(LST)
	counts 	= list()
	while j < N :
		if not single:
			start, stop = LST[j][0]-threshold, LST[j][1] + threshold
		else:
			start, stop = LST[j] -threshold, LST[j]  + threshold
			
		ct 	= 0
		if not single:
			while j < N and start < LST[j][1] and stop > LST[j][0]:
				start,stop 	= min(start, LST[j][0]-threshold), max(stop, LST[j][1]+threshold)
				ct+=1.0
				j+=1
		else:
			while j < N and start < LST[j]  and stop > LST[j] :
				start,stop 	= min(start, LST[j] -threshold), max(stop, LST[j] +threshold)
				ct+=1.0
				j+=1
			
		counts.append(ct)

		j+=1
	return counts
def draw_histogram(G,A):
	D_eRNA 	= [x for chrom in G for x in compute_distances(G[chrom])]
	D_ANNOT = [x for chrom in A for x in compute_distances(A[chrom])]
	D_NULL 	= [x for chrom in G for x in compute_distances_random(np.random.uniform(G[chrom][0][0], G[chrom][-1][1], int(len(G[chrom])))) ]
	F 		= plt.figure(figsize=(15,10))
	ax1 	= F.add_subplot(1,2,1)
	ax2 	= F.add_subplot(1,2,2)
	ax1.hist(D_eRNA, bins=75,alpha=0.3, normed=1, label="Bidir unannotated, N: " + str(sum([len(G[chrom]) for chrom in G])) )
	ax1.hist(D_NULL, bins=75, alpha=0.3, normed=1, label="NULL, N: "+ str(sum([len(G[chrom]) for chrom in G])))
	ax1.hist(D_ANNOT, bins=75, alpha=0.3, normed=1, label="TSS Bidir., N: "+ str(sum([len(A[chrom]) for chrom in A])))
	
	ax1.set_xlabel(r"$\log_{10}x_{t+1} - x_{t}$")
	ax1.set_ylabel("normalized frequency")
	ax1.grid()
	ax1.legend()
	thresholds 	= np.linspace(1000,10000000, 100)
	U 		= dict([(chrom,np.random.uniform(G[chrom][0][0], G[chrom][-1][1], int(len(G[chrom]) )))  for chrom in G])
	D_eRNA 	= [[x for chrom in G for x in linkage(G[chrom], threshold=t)  ]  for t in thresholds ]
	D_NULL 	= [[x for chrom in U for x in linkage(U[chrom], threshold=t, single=True)  ]  for t in thresholds ]
	ax2.scatter(thresholds, [np.mean(x) for x in D_eRNA], label="bidirectional", s=10)
	ax2.scatter(thresholds, [np.mean(x) for x in D_NULL], label="NULL", color="red", s=10)
	ax2.legend(loc=(0.6,0.5))
	ax2.set_xlabel("Linkage Threshold")
	ax2.set_ylabel("Number of Occupants")

	ax2.grid()

	plt.show()

if __name__ == "__main__":
	TSS_T 	= get_TSS_tree("/Users/joazofeifa/Lab/genome_files/TSS.bed")
	BIDIR 	= load("/Users/joazofeifa/Lab/gro_seq_files/Allen2014/EMG_out_files/Allen2014_boot_DMSO2_3-2_bootstrapped_bidirectional_hits_intervals.bed", bidir_bootstrap=True,FILTER=TSS_T)
	BIDIR_I = load("/Users/joazofeifa/Lab/gro_seq_files/Allen2014/EMG_out_files/Allen2014_boot_DMSO2_3-2_bootstrapped_bidirectional_hits_intervals.bed", bidir_bootstrap=True,isolate=TSS_T)
	
	draw_histogram(BIDIR, BIDIR_I)