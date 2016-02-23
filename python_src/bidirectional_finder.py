import numpy as np
import time
import model as MODEL
import matplotlib.pyplot as plt
import math
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
def window(X, win=10):
	i,j,k 				= 0,0,0
	scores, forward,reverse = list(), list(),list()
	while i < 5000:
		while j < X.shape[0] and  (X[j,0]-X[i,0]) < -win:
			j+=1
		while k <X.shape[0] and (X[k,0]-X[i,0]) < win:
			k+=1
		if i < X.shape[0] and j < X.shape[0] and k < X.shape[0]:
			a,b 					= X[j,0],X[k,0]
			N_forward,N_reverse  	= np.sum(X[j:k,1]),np.sum(X[j:k,2])
			N 						= N_forward + N_reverse
			pi 						= (N_forward+1) / (N+2)

			null 					= math.log(pi / (b-a))*N_forward +  math.log((1-pi) / (b-a))*N_reverse
			center 					= X[i,0]
			rvs 					= [MODEL.component_bidir(center, 1.0, 0.5, 0.9,pi , None) , MODEL.component_elongation( a,b, 0.1, pi, None, None, None, None, )]
			model 					= sum([math.log(sum([rv.pdf(X[u,0],1) for rv in rvs] ))*X[u,1] for u in range(j,k)])
			model 					+=sum([math.log(sum([rv.pdf(X[u,0],-1) for rv in rvs] ))*X[u,2] for u in range(j,k)])
			scores.append(model/null)


		i+=1
	return scores
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
		X[:,0]/=100.
		G[chrom]=X
	return G
def scan(G, win=1000):
	FHW 	= open("/Users/joazofeifa/test_peak_caller.bed","w")
	for chrom in G:
		print chrom, 
		scores 	= window(G[chrom]) 
		
		intervals  	= list()
		start 		= None
		for i in range( len(scores) ):
			if scores[i] < 1.0 and start is None:
				start 	= i
			elif start is not None and scores[i] > 1.0 :
				stop 	= i
				intervals.append((G[chrom][start,0]*100, G[chrom][stop,0]*100))
				start 	= None
		print len(intervals)
		for start, stop in intervals:
			FHW.write(chrom+"\t" + str(int(start)) + "\t" + str(int(stop)) + "\n")






if __name__== "__main__":
	forward 	= "/Users/joazofeifa/Lab/gro_seq_files/HCT116/bed_graph_files/test_DMSO2_3.pos.BedGraph"
	reverse 	= "/Users/joazofeifa/Lab/gro_seq_files/HCT116/bed_graph_files/test_DMSO2_3.neg.BedGraph"
	total,forward, reverse 		= get_counts((forward,reverse))
	G 			= bin(forward, reverse)
	scan(G)
	
