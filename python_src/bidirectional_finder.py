import numpy as np
import time
def get_counts(FILES):
	G= {}
	DS 	= [{}, {}]
	for i,FILE in enumerate(FILES):
		t=0
		with open(FILE) as FH:
			for line in FH:
				chrom,start, stop, cov 	= line.strip("\n").split("\t")
				if chrom not in DS[i]:
					t+=1
					if t > 1:
						break
					DS[i][chrom]=list()
				if chrom not in G:
					G[chrom] 	= [0,np.inf, -np.inf]
				x 				= (float(start) + float(stop))/2.
				G[chrom][0] 	+=abs(float(cov))
				G[chrom][1] 	= min((G[chrom][1], int(start)))
				G[chrom][2] 	= max((G[chrom][2], int(stop)))
				DS[i][chrom].append((x, float(cov)))
	return G,DS[0],DS[1]
def window(X, win=1000):
	i 				= 0
	forward,reverse = list(), list()
	while i < X.shape[0]:
		j,k 	= i,i
		while j < X.shape[0] and (X[j,0] - X[i,0]) < win:
			j+=1
		while k >=0 and (X[i,0] - X[k,0]) < win:
			k-=1
		if j < X.shape[0] and i < X.shape[0]:
			forward.append( (np.sum(X[i:j, 1])   ))
			reverse.append( (np.sum(X[k:i, 2])  ))
		i+=1
	return forward, reverse
def bin(Forward, Reverse,step_size=25):
	G 	= {}
	for chrom in Forward:
		forward, reverse 	= Forward[chrom], Reverse[chrom]
		start, stop 	= min((min(forward)[0],min(reverse)[0])),max((max(forward)[0],max(reverse)[0]))
		bins 			= int( (stop-start)/step_size)
		X 				= np.zeros((bins, 3))
		X[:,0] 			= np.linspace(start, stop, bins)
		for j,f in enumerate((forward, reverse)):
			i = 0
			for x,c in f:
				while i < X.shape[0] and X[i,0] <= x:
					i+=1
				X[i-1, j+1]+=c
		G[chrom]=X
	return G
def scan(G, win=1000):
	FHW 	= open("/Users/joazofeifa/test_peak_caller.bed","w")
	for chrom in G:
		windows_forward, windows_reverse 	= window(G[chrom], win=win) 
		Nf, Nr 		= sum( windows_forward), sum( windows_reverse)
		l 			= float(len(windows_forward)  )*win
		start 		= None
		scores 		= list()
		ef, er 		= Nf * (25.0 / l), Nr*(25.0 / l)
		
		intervals  	= list()
		for i in range( len(windows_forward) ):
			

			if windows_forward[i] > ef and windows_reverse[i] > er and start is None:
				start 	= i
			elif start is not None and (windows_forward[i] < ef or windows_reverse[i] < er) :
				stop 	= i
				intervals.append((start, stop))
				start 	= None

			elif (windows_forward[i] < ef or windows_reverse[i] < er) :
				start  	= None







if __name__== "__main__":
	forward 	= "/Users/joazofeifa/Lab/gro_seq_files/HCT116/bed_graph_files/test_DMSO2_3.pos.BedGraph"
	reverse 	= "/Users/joazofeifa/Lab/gro_seq_files/HCT116/bed_graph_files/test_DMSO2_3.neg.BedGraph"
	total,forward, reverse 		= get_counts((forward,reverse))
	G 			= bin(forward, reverse)
	scan(G)
	
